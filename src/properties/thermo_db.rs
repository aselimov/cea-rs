use crate::properties::{
    PropertiesError,
    error::make_parse_error,
    thermo_fit::{Phase, SpeciesElement, SpeciesThermoData, ThermoPolynomial},
    utils::parse_fields,
};

pub struct ThermoDB {
    pub products: Vec<SpeciesThermoData>,
    pub reactants: Vec<SpeciesThermoData>,
}

/// Parse a thermo formatted db
impl ThermoDB {
    pub fn parse(thermo_inp: &str) -> Result<Self, PropertiesError> {
        let mut lines = thermo_inp.lines();
        let mut products = Vec::new();
        let mut reactants = Vec::new();
        let mut parse_products = true;
        // Skip comments
        while let Some(line) = lines.next() {
            if line.trim().is_empty() || line.starts_with("!") {
                continue;
            } else if line.contains("thermo") {
                _ = lines.next().ok_or(PropertiesError::InvalidFile)?;
                continue;
            } else if line.contains("END PRODUCTS") {
                parse_products = false;
            } else if line.contains("END REACTANTS") {
                break;
            } else if parse_products {
                products.push(parse_species(line, &mut lines)?);
            } else {
                reactants.push(parse_species(line, &mut lines)?);
            }
        }

        Ok(ThermoDB {
            products,
            reactants,
        })
    }
}

fn parse_species<'a>(
    line: &str,
    lines: &mut impl Iterator<Item = &'a str>,
) -> Result<SpeciesThermoData, PropertiesError> {
    // Parsing a fortran generated file which means we used fixed column width parsing. Define the
    // fixed column widths used
    const SPECIES_LINE_2_WIDTHS: &[usize] = &[3, 7, 2, 6, 2, 6, 2, 6, 2, 6, 2, 6, 2, 13, 15];

    let name = line
        .get(0..16)
        .ok_or(PropertiesError::InvalidLine("name".to_string()))?
        .trim()
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

    // 0-interval species still have one reference state line (298.15 K data) that must be consumed
    if intervals == 0 {
        lines.next().ok_or(PropertiesError::InvalidFile)?;
    }

    Ok(SpeciesThermoData::new(
        name,
        elements,
        phase,
        polynomials,
        molecular_weight,
        h_formation,
    ))
}

fn parse_polynomials_block<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    intervals: usize,
) -> Result<Vec<ThermoPolynomial>, PropertiesError> {
    // Now parse the actual polynomial intervals
    (0..intervals)
        .map(|_| parse_polynomial_block(lines))
        .collect()
}

fn parse_polynomial_block<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
) -> Result<ThermoPolynomial, PropertiesError> {
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
        .map_err(|_| make_parse_error("temp_hi", "f64", &splits[1]))?;

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

    Ok(ThermoPolynomial {
        a,
        temp_range: (temp_lo, temp_hi),
    })
}

#[cfg(test)]
mod test {
    use crate::{
        assert_delta, assert_vec_delta,
        properties::{
            thermo_db::{ThermoDB, parse_polynomial_block, parse_polynomials_block, parse_species},
            thermo_fit::Phase,
        },
    };

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
        assert_delta!(polynomial.temp_range.1, 6000.000, 1e-3);
    }

    #[test]
    fn test_parse_polynomials_block() {
        let polynomials_block = r#"    300.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6918.671
 5.006608890D+03 1.861304407D+01 2.412531111D+00 1.987604647D-04-2.432362152D-07
 1.538281506D-10-3.944375734D-14                 3.887412680D+04 6.086585765D+00
   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6918.671
-2.920820938D+04 1.167751876D+02 2.356906505D+00 7.737231520D-05-1.529455262D-08
-9.971670260D-13 5.053278264D-16                 3.823288650D+04 6.600920155D+00
   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6918.671
-5.040682320D+08 3.802322650D+05-1.082347159D+02 1.549444292D-02-1.070103856D-06
 3.592110900D-11-4.696039394D-16                -2.901050501D+06 9.491883160D+02"#;

        let mut lines = polynomials_block.lines();
        let polynomials = parse_polynomials_block(&mut lines, 3).unwrap();

        let real_coeff_1 = [
            5.006608890e+03,
            1.861304407e+01,
            2.412531111e+00,
            1.987604647e-04,
            -2.432362152e-07,
            1.538281506e-10,
            -3.944375734e-14,
            3.887412680e+04,
            6.086585765e+00,
        ];

        assert_vec_delta!(real_coeff_1, polynomials[0].a, 1e-9);
        assert_delta!(polynomials[0].temp_range.0, 300.000, 1e-3);
        assert_delta!(polynomials[0].temp_range.1, 1000.000, 1e-3);

        let real_coeff_2 = [
            -2.920820938e+04,
            1.167751876e+02,
            2.356906505e+00,
            7.737231520e-05,
            -1.529455262e-08,
            -9.971670260e-13,
            5.053278264e-16,
            3.823288650e+04,
            6.600920155e+00,
        ];
        assert_vec_delta!(real_coeff_2, polynomials[1].a, 1e-9);
        assert_delta!(polynomials[1].temp_range.0, 1000.000, 1e-3);
        assert_delta!(polynomials[1].temp_range.1, 6000.000, 1e-3);

        let real_coeff_3 = [
            -5.040682320e+08,
            3.802322650e+05,
            -1.082347159e+02,
            1.549444292e-02,
            -1.070103856e-06,
            3.592110900e-11,
            -4.696039394e-16,
            -2.901050501e+06,
            9.491883160e+02,
        ];
        assert_vec_delta!(real_coeff_3, polynomials[2].a, 1e-9);
        assert_delta!(polynomials[2].temp_range.0, 6000.000, 1e-3);
        assert_delta!(polynomials[2].temp_range.1, 20000.000, 1e-3);
    }

    #[test]
    fn test_parse_species() {
        let species = r#"ALBr2             Gurvich,1996a pt1 p186 pt2 p149.
 2 tpis96 AL  1.00BR  2.00    0.00    0.00    0.00 0  186.7895380    -140662.125
    300.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        13397.875
 3.199375870D+04-7.119178970D+02 9.478258110D+00-4.875531670D-03 5.516512990D-06
-3.340053040D-09 8.368476840D-13                -1.540591306D+04-1.742171366D+01
   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        13397.875
-3.523782900D+05 4.671544170D+02 7.111908190D+00-5.551709200D-04 3.166301130D-07
-5.521028330D-11 3.176725950D-15                -2.265004078D+04-2.695610360D+00"#;

        let mut lines = species.lines();
        let line = lines.next().unwrap();
        let species = parse_species(line, &mut lines).unwrap();

        assert_eq!(species.name, "ALBr2");
        assert_eq!(species.elements.len(), 2);
        assert_eq!(species.elements[0].element, "AL");
        assert_eq!(species.elements[0].count, 1.0);
        assert_eq!(species.elements[1].element, "BR");
        assert_eq!(species.elements[1].count, 2.0);
        assert!(matches!(species.phase, Phase::Gas));
        assert_delta!(species.molecular_weight, 186.7895380, 1e-7);
        assert_delta!(species.h_formation, -140662.125, 1e-3);

        let real_coeff_1 = [
            3.199375870e+04,
            -7.119178970e+02,
            9.478258110e+00,
            -4.875531670e-03,
            5.516512990e-06,
            -3.340053040e-09,
            8.368476840e-13,
            -1.540591306e+04,
            -1.742171366e+01,
        ];

        assert_vec_delta!(species.polynomial_at(650.0).a, real_coeff_1, 1e-9);
        assert_delta!(species.polynomial_at(650.0).temp_range.0, 300.000, 1e-3);
        assert_delta!(species.polynomial_at(650.0).temp_range.1, 1000.000, 1e-3);

        let real_coeff_2 = [
            -3.523782900e+05,
            4.671544170e+02,
            7.111908190e+00,
            -5.551709200e-04,
            3.166301130e-07,
            -5.521028330e-11,
            3.176725950e-15,
            -2.265004078e+04,
            -2.695610360e+00,
        ];

        assert_vec_delta!(species.polynomial_at(3500.0).a, real_coeff_2, 1e-9);
        assert_delta!(species.polynomial_at(3500.0).temp_range.0, 1000.000, 1e-3);
        assert_delta!(species.polynomial_at(3500.0).temp_range.1, 6000.000, 1e-3);
    }

    #[test]
    fn test_parse_thermo_db() {
        let thermo_file_contents = r#"!
! Some pointless header lines
!

thermo
    200.00   1000.00   6000.00  20000.   9/8/2021
ALCL3             Gurvich,1996a pt1 p173 pt2 p134.
 2 tpis96 AL  1.00CL  3.00    0.00    0.00    0.00 0  133.3405380    -584678.863
    300.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        16400.803
 7.750600970D+04-1.440779717D+03 1.401744141D+01-6.381631240D-03 5.871674720D-06
-2.908872278D-09 5.994050890D-13                -6.579343180D+04-4.494017799D+01
   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        16400.803
-1.378630916D+05-5.579207290D+01 1.004190387D+01-1.682165339D-05 3.724664660D-09
-4.275526780D-13 1.982341329D-17                -7.343407470D+04-2.045130429D+01
END PRODUCTS

Air               Mole%:N2 78.084,O2 20.9476,Ar .9365,CO2 .0319.Gordon,1982.Reac
 2 g 9/95 N 1.5617O .41959AR.00937C .00032  .00000 0   28.9651159       -125.530
    300.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8649.264
 1.009950160D+04-1.968275610D+02 5.009155110D+00-5.761013730D-03 1.066859930D-05
-7.940297970D-09 2.185231910D-12                -1.767967310D+02-3.921504225D+00
   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8649.264
 2.415214430D+05-1.257874600D+03 5.144558670D+00-2.138541790D-04 7.065227840D-08
-1.071483490D-11 6.577800150D-16                 6.462263190D+03-8.147411905D+00
n-Butanol         ANL's Active Thermochemical Tables (ATcT).              React.
 0 g 5/23 C   4.00H  10.00O   1.00   0.00     0.00 1   74.1216000    -278510.000
    298.150      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000
END REACTANTS
"#;
        let thermo_db = ThermoDB::parse(thermo_file_contents).unwrap();

        assert_eq!(thermo_db.products.len(), 1);
        assert_eq!(thermo_db.reactants.len(), 2);

        // --- ALCL3 (product) ---
        let alcl3 = &thermo_db.products[0];
        assert_eq!(alcl3.name, "ALCL3");
        assert_eq!(alcl3.elements.len(), 2);
        assert_eq!(alcl3.elements[0].element, "AL");
        assert_delta!(alcl3.elements[0].count, 1.0, 1e-9);
        assert_eq!(alcl3.elements[1].element, "CL");
        assert_delta!(alcl3.elements[1].count, 3.0, 1e-9);
        assert!(matches!(alcl3.phase, Phase::Gas));
        assert_delta!(alcl3.molecular_weight, 133.3405380, 1e-7);
        assert_delta!(alcl3.h_formation, -584678.863, 1e-3);
        assert_eq!(alcl3.num_polynomials(), 2);

        assert_vec_delta!(
            alcl3.polynomial_at(650.0).a,
            [
                7.750600970e+04,
                -1.440779717e+03,
                1.401744141e+01,
                -6.381631240e-03,
                5.871674720e-06,
                -2.908872278e-09,
                5.994050890e-13,
                -6.579343180e+04,
                -4.494017799e+01,
            ],
            1e-9
        );
        assert_delta!(alcl3.polynomial_at(650.0).temp_range.0, 300.0, 1e-3);
        assert_delta!(alcl3.polynomial_at(650.0).temp_range.1, 1000.0, 1e-3);

        assert_vec_delta!(
            alcl3.polynomial_at(3500.0).a,
            [
                -1.378630916e+05,
                -5.579207290e+01,
                1.004190387e+01,
                -1.682165339e-05,
                3.724664660e-09,
                -4.275526780e-13,
                1.982341329e-17,
                -7.343407470e+04,
                -2.045130429e+01,
            ],
            1e-9
        );
        assert_delta!(alcl3.polynomial_at(3500.0).temp_range.0, 1000.0, 1e-3);
        assert_delta!(alcl3.polynomial_at(3500.0).temp_range.1, 6000.0, 1e-3);

        // --- Air (reactant 0) ---
        let air = &thermo_db.reactants[0];
        assert_eq!(air.name, "Air");
        assert_eq!(air.elements.len(), 4);
        assert_eq!(air.elements[0].element, "N");
        assert_delta!(air.elements[0].count, 1.5617, 1e-9);
        assert_eq!(air.elements[1].element, "O");
        assert_delta!(air.elements[1].count, 0.41959, 1e-9);
        assert_eq!(air.elements[2].element, "AR");
        assert_delta!(air.elements[2].count, 0.00937, 1e-9);
        assert_eq!(air.elements[3].element, "C");
        assert_delta!(air.elements[3].count, 0.00032, 1e-9);
        assert!(matches!(air.phase, Phase::Gas));
        assert_delta!(air.molecular_weight, 28.9651159, 1e-7);
        assert_delta!(air.h_formation, -125.530, 1e-3);
        assert_eq!(air.num_polynomials(), 2);

        assert_vec_delta!(
            air.polynomial_at(650.0).a,
            [
                1.009950160e+04,
                -1.968275610e+02,
                5.009155110e+00,
                -5.761013730e-03,
                1.066859930e-05,
                -7.940297970e-09,
                2.185231910e-12,
                -1.767967310e+02,
                -3.921504225e+00,
            ],
            1e-9
        );
        assert_delta!(air.polynomial_at(650.0).temp_range.0, 300.0, 1e-3);
        assert_delta!(air.polynomial_at(650.0).temp_range.1, 1000.0, 1e-3);

        assert_vec_delta!(
            air.polynomial_at(3500.0).a,
            [
                2.415214430e+05,
                -1.257874600e+03,
                5.144558670e+00,
                -2.138541790e-04,
                7.065227840e-08,
                -1.071483490e-11,
                6.577800150e-16,
                6.462263190e+03,
                -8.147411905e+00,
            ],
            1e-9
        );
        assert_delta!(air.polynomial_at(3500.0).temp_range.0, 1000.0, 1e-3);
        assert_delta!(air.polynomial_at(3500.0).temp_range.1, 6000.0, 1e-3);

        // --- n-Butanol (reactant 1) ---
        let butanol = &thermo_db.reactants[1];
        assert_eq!(butanol.name, "n-Butanol");
        assert_eq!(butanol.elements.len(), 3);
        assert_eq!(butanol.elements[0].element, "C");
        assert_delta!(butanol.elements[0].count, 4.0, 1e-9);
        assert_eq!(butanol.elements[1].element, "H");
        assert_delta!(butanol.elements[1].count, 10.0, 1e-9);
        assert_eq!(butanol.elements[2].element, "O");
        assert_delta!(butanol.elements[2].count, 1.0, 1e-9);
        assert!(matches!(butanol.phase, Phase::Condensed));
        assert_delta!(butanol.molecular_weight, 74.1216000, 1e-7);
        assert_delta!(butanol.h_formation, -278510.000, 1e-3);
        assert_eq!(butanol.num_polynomials(), 0);
    }
}
