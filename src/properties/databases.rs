use crate::properties::{
    PropertiesError,
    error::make_parse_error,
    polynomials::{Phase, Polynomial, SpeciesElement, SpeciesPolynomial},
};

fn parse_fields<'a>(line: &'a str, widths: &[usize]) -> Vec<&'a str> {
    let mut fields = Vec::new();
    let mut pos = 0;

    for &width in widths {
        if let Some(field) = line.get(pos..pos + width) {
            fields.push(field.trim()); // .trim() removes padding spaces
        }
        pos += width;
    }

    fields
}

/// Parse a thermo formatted db
fn parse_thermo(thermo_inp: &str) -> Result<(), PropertiesError> {
    // Parsing a fortran generated file which means we used fixed column width parsing. Define the
    // fixed column widths used
    const SPECIES_LINE_2_WIDTHS: &[usize] = &[3, 7, 2, 6, 2, 6, 2, 6, 2, 6, 2, 6, 2, 13, 15];
    // Ignore the coefficients since they are the same
    const SPECIES_INTERVAL_1_WIDTHS: &[usize] = &[11, 11];
    const SPECIES_INTERVAL_2_WIDTHS: &[usize] = &[16; 5];
    const SPECIES_INTERVAL_3_WIDTHS: &[usize] = &[16; 5];
    let mut lines = thermo_inp.lines();
    let mut name;
    let mut species = Vec::new();
    let mut parse_species = true;
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
        if parse_species {
            //Line 1 name + metadata
            name = line
                .get(0..16)
                .ok_or(PropertiesError::InvalidLine("name".to_string()))?
                .to_string();

            // line 2
            let split = parse_fields(line, SPECIES_LINE_2_WIDTHS);
            let intervals: usize = split[0]
                .parse()
                .map_err(|_| make_parse_error("intervals", "usize", split[0]))?;

            let mut elements = vec![];
            for i in (2..=10).step_by(2) {
                let element = split[i].to_string();
                let count: f64 = split[i + 1]
                    .parse()
                    .map_err(|_| make_parse_error("species_count", "f64", split[i + 1]))?;

                if count.abs() > 1e-8 {
                    elements.push(SpeciesElement { element, count })
                }
            }

            let phase = match split[12]
                .parse::<i32>()
                .map_err(|_| make_parse_error("phase", "i32", split[12]))?
            {
                0 => Phase::Gas,
                _ => Phase::Condensed,
            };

            let molecular_weight = split[13]
                .parse()
                .map_err(|_| make_parse_error("molecular_weight", "f64", split[13]))?;

            let h_formation = split[14]
                .parse()
                .map_err(|_| make_parse_error("h_formation", "f64", split[14]))?;

            // Now parse the actual polynomial intervals
            let polynomials = (0..intervals)
                .map(|_| {
                    // Parse only the temps from first line
                    let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
                    let splits = parse_fields(line, SPECIES_INTERVAL_1_WIDTHS);
                    let temp_lo: f64 = splits[0]
                        .parse()
                        .map_err(|_| make_parse_error("temp_lo", "f64", splits[0]))?;

                    let temp_hi: f64 = splits[1]
                        .parse()
                        .map_err(|_| make_parse_error("temp_hi", "f64", splits[0]))?;

                    // Now parse the first 5 coefficients
                    let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
                    let splits = parse_fields(line, SPECIES_INTERVAL_2_WIDTHS);
                    let mut a: Vec<f64> = splits
                        .iter()
                        .map(|val| val.parse().map_err(|_| make_parse_error("a", "f64", val)))
                        .collect::<Result<Vec<f64>, PropertiesError>>()?;

                    let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
                    let splits = parse_fields(line, SPECIES_INTERVAL_3_WIDTHS);
                    //TODO: Finish parsing a. Just need to parse the 4 values without the empty space.
                    //It's o6 a7 <space> b1 b2
                    for i in [0, 1, 3, 4] {
                        a.push(
                            splits[i]
                                .parse()
                                .map_err(|_| make_parse_error("a", "f64", splits[i]))?,
                        );
                    }

                    Ok(Polynomial {
                        a,
                        temp_range: (temp_lo, temp_hi),
                    })
                })
                .collect::<Result<Vec<Polynomial>, PropertiesError>>()?;

            species.push(SpeciesPolynomial {
                name,
                polynomials,
                elements,
                phase,
                molecular_weight,
                h_formation,
            });
        }
    }
    todo!()
}
