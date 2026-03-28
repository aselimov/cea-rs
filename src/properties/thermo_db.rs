use crate::properties::{
    PropertiesError,
    error::make_parse_error,
    polynomials::{Phase, Polynomial, SpeciesElement, SpeciesPolynomial},
};

fn parse_fields<'a>(line: &'a str, widths: &[usize]) -> Vec<String> {
    let mut fields = Vec::new();
    let mut pos = 0;

    for &width in widths {
        if let Some(field) = line.get(pos..pos + width) {
            // The replace chnages the fortran formatted D exponential for the normal E exponential
            fields.push(field.trim().replace("D", "E"));
        }
        pos += width;
    }

    fields
}

pub struct ThermoDB {
    pub products: Vec<SpeciesPolynomial>,
    pub reactants: Vec<SpeciesPolynomial>,
}

/// Parse a thermo formatted db
impl ThermoDB {
    pub fn parse(thermo_inp: &str) -> Result<Self, PropertiesError> {
        let mut lines = thermo_inp.lines();
        let mut species = Vec::new();
        let species_block = true;
        // Skip comments
        while let Some(line) = lines.next() {
            if line.starts_with("!") {
                continue;
            }

            // Skip pointless header lines
            if line.contains("thermo") {
                _ = lines.next().ok_or(PropertiesError::InvalidFile)?;
                continue;
            }

            // Parse species block
            if species_block {
                species.push(Self::parse_species(&mut lines)?);
            }
            //TODO Distinguish between products and reactants
        }
        todo!()
    }

    fn parse_species<'a>(
        lines: &mut impl Iterator<Item = &'a str>,
    ) -> Result<SpeciesPolynomial, PropertiesError> {
        // Parsing a fortran generated file which means we used fixed column width parsing. Define the
        // fixed column widths used
        const SPECIES_LINE_2_WIDTHS: &[usize] = &[3, 7, 2, 6, 2, 6, 2, 6, 2, 6, 2, 6, 2, 13, 15];

        let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
        let name = line
            .get(0..16)
            .ok_or(PropertiesError::InvalidLine("name".to_string()))?
            .to_string();

        // line 2
        let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
        let split = parse_fields(line, SPECIES_LINE_2_WIDTHS);
        let intervals: usize = split[0]
            .parse()
            .map_err(|_| make_parse_error("intervals", "usize", &split[0]))?;

        let mut elements = vec![];
        for i in (2..=10).step_by(2) {
            let element = split[i].to_string();
            let count: f64 = split[i + 1]
                .parse()
                .map_err(|_| make_parse_error("species_count", "f64", &split[i + 1]))?;

            if count.abs() > 1e-8 {
                elements.push(SpeciesElement { element, count })
            }
        }

        let phase = match split[12]
            .parse::<i32>()
            .map_err(|_| make_parse_error("phase", "i32", &split[12]))?
        {
            0 => Phase::Gas,
            _ => Phase::Condensed,
        };

        let molecular_weight = split[13]
            .parse()
            .map_err(|_| make_parse_error("molecular_weight", "f64", &split[13]))?;

        let h_formation = split[14]
            .parse()
            .map_err(|_| make_parse_error("h_formation", "f64", &split[14]))?;

        let polynomials = parse_polynomials_block(lines, intervals)?;

        Ok(SpeciesPolynomial {
            name,
            polynomials,
            elements,
            phase,
            molecular_weight,
            h_formation,
        })
    }
}

fn parse_polynomials_block<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    intervals: usize,
) -> Result<Vec<Polynomial>, PropertiesError> {
    // Now parse the actual polynomial intervals
    (0..intervals)
        .map(|_| parse_polynomial_block(lines))
        .collect()
}

fn parse_polynomial_block<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
) -> Result<Polynomial, PropertiesError> {
    // Ignore the coefficients since they are the same
    const SPECIES_INTERVAL_1_WIDTHS: &[usize] = &[11, 11];
    const SPECIES_INTERVAL_2_WIDTHS: &[usize] = &[16; 5];
    const SPECIES_INTERVAL_3_WIDTHS: &[usize] = &[16; 5];

    // Parse only the temps from first line
    let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
    let splits = parse_fields(line, SPECIES_INTERVAL_1_WIDTHS);
    let temp_lo: f64 = splits[0]
        .parse()
        .map_err(|_| make_parse_error("temp_lo", "f64", &splits[0]))?;

    let temp_hi: f64 = splits[1]
        .parse()
        .map_err(|_| make_parse_error("temp_hi", "f64", &splits[0]))?;

    // Now parse the first 5 coefficients
    let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
    let splits = parse_fields(line, SPECIES_INTERVAL_2_WIDTHS);
    let mut a: Vec<f64> = splits
        .iter()
        .map(|val| val.parse().map_err(|_| make_parse_error("a", "f64", val)))
        .collect::<Result<Vec<f64>, PropertiesError>>()?;

    let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
    let splits = parse_fields(line, SPECIES_INTERVAL_3_WIDTHS);
    for i in [0, 1, 3, 4] {
        a.push(
            splits[i]
                .parse()
                .map_err(|_| make_parse_error("a", "f64", &splits[i]))?,
        );
    }

    Ok(Polynomial {
        a,
        temp_range: (temp_lo, temp_hi),
    })
}

#[cfg(test)]
mod test {
    use crate::{assert_delta, assert_vec_delta, properties::thermo_db::parse_polynomial_block};

    #[test]
    fn test_parse_polynomial_block() {
        let polynomial_block = r#"   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428
-4.181183250D+03-9.948557270D+00 2.548615878D+00-5.878760040D-05 3.132291294D-08
-7.748894630D-12 7.274447690D-16                 1.091011485D+05 3.488667290D+00"#;
        let mut lines = polynomial_block.lines();
        let polynomial = parse_polynomial_block(&mut lines).unwrap();

        let real = [
            -4.181183250e+03,
            -9.948557270e+00,
            2.548615878e+00,
            -5.878760040e-05,
            3.132291294e-08,
            -7.748894630e-12,
            7.274447690e-16,
            1.091011485e+05,
            3.488667290e+00,
        ];

        assert_vec_delta!(real, polynomial.a, 1e-9);
        assert_delta!(polynomial.temp_range.0, 1000.000, 1e-3);
        assert_delta!(polynomial.temp_range.1, 6000.0007, 1e-3);
    }
}
