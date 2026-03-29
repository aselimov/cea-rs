use crate::properties::{
    PropertiesError,
    error::make_parse_error,
    transport_fit::{SpeciesTransportData, TransportFit},
    utils::parse_fields,
};

enum ViscosityOrConductivity {
    Viscosity,
    Conductivity,
}

pub struct TransportDB {
    data: Vec<SpeciesTransportData>,
}

impl TransportDB {
    pub fn parse(transport_file_contents: &str) -> Result<Self, PropertiesError> {
        let mut lines = transport_file_contents.lines().peekable();
        let mut data = vec![];
        // Ignore header line
        _ = lines.next();
        loop {
            data.push(parse_species_transport_block(&mut lines)?);

            if let Some(line) = lines.peek() {
                if line.contains("end") {
                    break;
                }
            } else {
                break;
            }
        }

        Ok(TransportDB { data })
    }
}

fn parse_species_transport_block<'a>(
    mut lines: impl Iterator<Item = &'a str>,
) -> Result<SpeciesTransportData, PropertiesError> {
    let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
    let mut viscosities = vec![];
    let mut conductivities = vec![];
    let (name, viscosity_count, conductivity_count) = parse_species_header_line(line)?;
    for _ in 0..viscosity_count + conductivity_count {
        let line = lines.next().ok_or(PropertiesError::InvalidFile)?;
        match parse_coefficients(line)? {
            (ViscosityOrConductivity::Viscosity, fit) => viscosities.push(fit),
            (ViscosityOrConductivity::Conductivity, fit) => conductivities.push(fit),
        }
    }
    Ok(SpeciesTransportData {
        name,
        viscosities,
        conductivities,
    })
}

fn parse_species_header_line(line: &str) -> Result<(String, usize, usize), PropertiesError> {
    const HEADER_WIDTHS: &[usize] = &[34, 1, 1, 1, 1];

    let splits = parse_fields(line, HEADER_WIDTHS);

    let name = splits[0].to_string();
    let viscosity_count: usize = splits[2]
        .parse()
        .map_err(|_| make_parse_error("V", "usize", &splits[2]))?;

    let conductivity_count: usize = splits[4]
        .parse()
        .map_err(|_| make_parse_error("C", "usize", &splits[2]))?;
    Ok((name, viscosity_count, conductivity_count))
}

fn parse_coefficients(
    line: &str,
) -> Result<(ViscosityOrConductivity, TransportFit), PropertiesError> {
    const COEFFICIENTS_WIDTH: &[usize] = &[2, 7, 11, 15, 15, 15, 15];
    let splits = parse_fields(line, COEFFICIENTS_WIDTH);
    let transport_type = if splits[0].contains("V") {
        ViscosityOrConductivity::Viscosity
    } else {
        ViscosityOrConductivity::Conductivity
    };
    Ok((
        transport_type,
        TransportFit {
            temp_range: (
                splits[1]
                    .parse()
                    .map_err(|_| make_parse_error("tlo", "f64", &splits[1]))?,
                splits[2]
                    .parse()
                    .map_err(|_| make_parse_error("tlo", "f64", &splits[2]))?,
            ),
            a: splits[3]
                .parse()
                .map_err(|_| make_parse_error("a", "f64", &splits[3]))?,
            b: splits[4]
                .parse()
                .map_err(|_| make_parse_error("a", "f64", &splits[4]))?,
            c: splits[5]
                .parse()
                .map_err(|_| make_parse_error("a", "f64", &splits[5]))?,

            d: splits[6]
                .parse()
                .map_err(|_| make_parse_error("a", "f64", &splits[6]))?,
        },
    ))
}
#[cfg(test)]
pub mod test {
    use crate::{
        assert_delta,
        properties::transport_db::{
            TransportDB, ViscosityOrConductivity, parse_coefficients, parse_species_header_line,
            parse_species_transport_block,
        },
    };

    #[test]
    fn test_parse_species_transport_header() {
        let line = "Ar                                V3C1  BICH ET AL (1990)";

        let (name, viscosity_count, conductivity_count) = parse_species_header_line(line).unwrap();

        assert_eq!(name, "Ar");
        assert_eq!(viscosity_count, 3);
        assert_eq!(conductivity_count, 1);
    }

    #[test]
    fn test_parse_coefficients() {
        let line =
            " V  200.0   1000.0   0.61205763E 00-0.67714354E 02 0.19040660E 03 0.21588272E 01";
        let (transport_type, fit) = parse_coefficients(line).unwrap();

        assert!(matches!(transport_type, ViscosityOrConductivity::Viscosity));
        assert_delta!(fit.temp_range.0, 200.0, 1e-1);
        assert_delta!(fit.temp_range.1, 1000.0, 1e-1);
        assert_delta!(fit.a, 0.61205763e00, 1e-8);
        assert_delta!(fit.b, -0.67714354E02, 1e-8);
        assert_delta!(fit.c, 0.19040660e03, 1e-8);
        assert_delta!(fit.d, 0.21588272E01, 1e-8);

        let line =
            " C 5000.0  15000.0   0.76269502E+00 0.62341752E+03-0.71899552E+06 0.56927918E+00";

        let (transport_type, fit) = parse_coefficients(line).unwrap();
        assert!(matches!(
            transport_type,
            ViscosityOrConductivity::Conductivity
        ));
        assert_delta!(fit.temp_range.0, 5000.0, 1e-1);
        assert_delta!(fit.temp_range.1, 15000.0, 1e-1);
        assert_delta!(fit.a, 0.76269502E+00, 1e-8);
        assert_delta!(fit.b, 0.62341752E+03, 1e-8);
        assert_delta!(fit.c, -0.71899552E+06, 1e-8);
        assert_delta!(fit.d, 0.56927918E+00, 1e-8);
    }

    #[test]
    fn test_parse_transport_species_block() {
        let species_block = r#"Ar                                V3C3  BICH ET AL (1990)
 V  200.0   1000.0   0.61205763E 00-0.67714354E 02 0.19040660E 03 0.21588272E 01
 V 1000.0   5000.0   0.69357334E 00 0.70953943E 02-0.28386007E 05 0.14856447E 01
 V 5000.0  15000.0   0.76608935E+00 0.67867215E+03-0.84991417E+06 0.77935167E+00
 C  200.0   1000.0   0.60968928E 00-0.70892249E 02 0.58420624E 03 0.19337152E 01
 C 1000.0   5000.0   0.69075463E 00 0.62676058E 02-0.25667413E 05 0.12664189E 01
 C 5000.0  15000.0   0.76269502E+00 0.62341752E+03-0.71899552E+06 0.56927918E+00"#;

        let lines = species_block.lines();
        let species_data = parse_species_transport_block(lines).unwrap();

        assert_eq!(species_data.name, "Ar");
        assert_eq!(species_data.viscosities.len(), 3);
        assert_eq!(species_data.conductivities.len(), 3);

        // Viscosity fits
        assert_delta!(species_data.viscosities[0].temp_range.0, 200.0, 1e-1);
        assert_delta!(species_data.viscosities[0].temp_range.1, 1000.0, 1e-1);
        assert_delta!(species_data.viscosities[0].a, 0.61205763e00, 1e-8);
        assert_delta!(species_data.viscosities[0].b, -0.67714354e02, 1e-8);
        assert_delta!(species_data.viscosities[0].c, 0.19040660e03, 1e-8);
        assert_delta!(species_data.viscosities[0].d, 0.21588272e01, 1e-8);

        assert_delta!(species_data.viscosities[1].temp_range.0, 1000.0, 1e-1);
        assert_delta!(species_data.viscosities[1].temp_range.1, 5000.0, 1e-1);
        assert_delta!(species_data.viscosities[1].a, 0.69357334e00, 1e-8);
        assert_delta!(species_data.viscosities[1].b, 0.70953943e02, 1e-8);
        assert_delta!(species_data.viscosities[1].c, -0.28386007e05, 1e-8);
        assert_delta!(species_data.viscosities[1].d, 0.14856447e01, 1e-8);

        assert_delta!(species_data.viscosities[2].temp_range.0, 5000.0, 1e-1);
        assert_delta!(species_data.viscosities[2].temp_range.1, 15000.0, 1e-1);
        assert_delta!(species_data.viscosities[2].a, 0.76608935e00, 1e-8);
        assert_delta!(species_data.viscosities[2].b, 0.67867215e03, 1e-8);
        assert_delta!(species_data.viscosities[2].c, -0.84991417e06, 1e-8);
        assert_delta!(species_data.viscosities[2].d, 0.77935167e00, 1e-8);

        // Conductivity fits
        assert_delta!(species_data.conductivities[0].temp_range.0, 200.0, 1e-1);
        assert_delta!(species_data.conductivities[0].temp_range.1, 1000.0, 1e-1);
        assert_delta!(species_data.conductivities[0].a, 0.60968928e00, 1e-8);
        assert_delta!(species_data.conductivities[0].b, -0.70892249e02, 1e-8);
        assert_delta!(species_data.conductivities[0].c, 0.58420624e03, 1e-8);
        assert_delta!(species_data.conductivities[0].d, 0.19337152e01, 1e-8);

        assert_delta!(species_data.conductivities[1].temp_range.0, 1000.0, 1e-1);
        assert_delta!(species_data.conductivities[1].temp_range.1, 5000.0, 1e-1);
        assert_delta!(species_data.conductivities[1].a, 0.69075463e00, 1e-8);
        assert_delta!(species_data.conductivities[1].b, 0.62676058e02, 1e-8);
        assert_delta!(species_data.conductivities[1].c, -0.25667413e05, 1e-8);
        assert_delta!(species_data.conductivities[1].d, 0.12664189e01, 1e-8);

        assert_delta!(species_data.conductivities[2].temp_range.0, 5000.0, 1e-1);
        assert_delta!(species_data.conductivities[2].temp_range.1, 15000.0, 1e-1);
        assert_delta!(species_data.conductivities[2].a, 0.76269502e00, 1e-8);
        assert_delta!(species_data.conductivities[2].b, 0.62341752e03, 1e-8);
        assert_delta!(species_data.conductivities[2].c, -0.71899552e06, 1e-8);
        assert_delta!(species_data.conductivities[2].d, 0.56927918e00, 1e-8);
    }

    #[test]
    fn test_parse_transport_file() {
        let trans_input = r#"transport property coefficients         
Ar                                V3C3  BICH ET AL (1990)
 V  200.0   1000.0   0.61205763E 00-0.67714354E 02 0.19040660E 03 0.21588272E 01
 V 1000.0   5000.0   0.69357334E 00 0.70953943E 02-0.28386007E 05 0.14856447E 01
 V 5000.0  15000.0   0.76608935E+00 0.67867215E+03-0.84991417E+06 0.77935167E+00
 C  200.0   1000.0   0.60968928E 00-0.70892249E 02 0.58420624E 03 0.19337152E 01
 C 1000.0   5000.0   0.69075463E 00 0.62676058E 02-0.25667413E 05 0.12664189E 01
 C 5000.0  15000.0   0.76269502E+00 0.62341752E+03-0.71899552E+06 0.56927918E+00
BCL3                              V2C2  SVEHLA (1962)
 V  300.0   1000.0   0.52572590E 00-0.27803504E 03 0.19159256E 05 0.24373790E 01
 V 1000.0   5000.0   0.62929553E 00-0.60723560E 02-0.37711618E 05 0.15615047E 01
 C  300.0   1000.0   0.41518585E 00-0.48149960E 03 0.30788060E 05 0.33168239E 01
 C 1000.0   5000.0   0.61148589E 00-0.18167042E 03-0.20976969E 05 0.17127671E 01
 end"#;

        let db = TransportDB::parse(trans_input).unwrap();

        assert_eq!(db.data.len(), 2);

        // Ar species (V3C3)
        let ar = &db.data[0];
        assert_eq!(ar.name, "Ar");
        assert_eq!(ar.viscosities.len(), 3);
        assert_eq!(ar.conductivities.len(), 3);

        assert_delta!(ar.viscosities[0].temp_range.0, 200.0, 1e-1);
        assert_delta!(ar.viscosities[0].temp_range.1, 1000.0, 1e-1);
        assert_delta!(ar.viscosities[0].a, 0.61205763e00, 1e-8);
        assert_delta!(ar.viscosities[0].b, -0.67714354e02, 1e-8);
        assert_delta!(ar.viscosities[0].c, 0.19040660e03, 1e-8);
        assert_delta!(ar.viscosities[0].d, 0.21588272e01, 1e-8);

        assert_delta!(ar.viscosities[1].temp_range.0, 1000.0, 1e-1);
        assert_delta!(ar.viscosities[1].temp_range.1, 5000.0, 1e-1);
        assert_delta!(ar.viscosities[1].a, 0.69357334e00, 1e-8);
        assert_delta!(ar.viscosities[1].b, 0.70953943e02, 1e-8);
        assert_delta!(ar.viscosities[1].c, -0.28386007e05, 1e-8);
        assert_delta!(ar.viscosities[1].d, 0.14856447e01, 1e-8);

        assert_delta!(ar.viscosities[2].temp_range.0, 5000.0, 1e-1);
        assert_delta!(ar.viscosities[2].temp_range.1, 15000.0, 1e-1);
        assert_delta!(ar.viscosities[2].a, 0.76608935e00, 1e-8);
        assert_delta!(ar.viscosities[2].b, 0.67867215e03, 1e-8);
        assert_delta!(ar.viscosities[2].c, -0.84991417e06, 1e-8);
        assert_delta!(ar.viscosities[2].d, 0.77935167e00, 1e-8);

        assert_delta!(ar.conductivities[0].temp_range.0, 200.0, 1e-1);
        assert_delta!(ar.conductivities[0].temp_range.1, 1000.0, 1e-1);
        assert_delta!(ar.conductivities[0].a, 0.60968928e00, 1e-8);
        assert_delta!(ar.conductivities[0].b, -0.70892249e02, 1e-8);
        assert_delta!(ar.conductivities[0].c, 0.58420624e03, 1e-8);
        assert_delta!(ar.conductivities[0].d, 0.19337152e01, 1e-8);

        assert_delta!(ar.conductivities[1].temp_range.0, 1000.0, 1e-1);
        assert_delta!(ar.conductivities[1].temp_range.1, 5000.0, 1e-1);
        assert_delta!(ar.conductivities[1].a, 0.69075463e00, 1e-8);
        assert_delta!(ar.conductivities[1].b, 0.62676058e02, 1e-8);
        assert_delta!(ar.conductivities[1].c, -0.25667413e05, 1e-8);
        assert_delta!(ar.conductivities[1].d, 0.12664189e01, 1e-8);

        assert_delta!(ar.conductivities[2].temp_range.0, 5000.0, 1e-1);
        assert_delta!(ar.conductivities[2].temp_range.1, 15000.0, 1e-1);
        assert_delta!(ar.conductivities[2].a, 0.76269502e00, 1e-8);
        assert_delta!(ar.conductivities[2].b, 0.62341752e03, 1e-8);
        assert_delta!(ar.conductivities[2].c, -0.71899552e06, 1e-8);
        assert_delta!(ar.conductivities[2].d, 0.56927918e00, 1e-8);

        // BCL3 species (V2C2)
        let bcl3 = &db.data[1];
        assert_eq!(bcl3.name, "BCL3");
        assert_eq!(bcl3.viscosities.len(), 2);
        assert_eq!(bcl3.conductivities.len(), 2);

        assert_delta!(bcl3.viscosities[0].temp_range.0, 300.0, 1e-1);
        assert_delta!(bcl3.viscosities[0].temp_range.1, 1000.0, 1e-1);
        assert_delta!(bcl3.viscosities[0].a, 0.52572590e00, 1e-8);
        assert_delta!(bcl3.viscosities[0].b, -0.27803504e03, 1e-8);
        assert_delta!(bcl3.viscosities[0].c, 0.19159256e05, 1e-8);
        assert_delta!(bcl3.viscosities[0].d, 0.24373790e01, 1e-8);

        assert_delta!(bcl3.viscosities[1].temp_range.0, 1000.0, 1e-1);
        assert_delta!(bcl3.viscosities[1].temp_range.1, 5000.0, 1e-1);
        assert_delta!(bcl3.viscosities[1].a, 0.62929553e00, 1e-8);
        assert_delta!(bcl3.viscosities[1].b, -0.60723560e02, 1e-8);
        assert_delta!(bcl3.viscosities[1].c, -0.37711618e05, 1e-8);
        assert_delta!(bcl3.viscosities[1].d, 0.15615047e01, 1e-8);

        assert_delta!(bcl3.conductivities[0].temp_range.0, 300.0, 1e-1);
        assert_delta!(bcl3.conductivities[0].temp_range.1, 1000.0, 1e-1);
        assert_delta!(bcl3.conductivities[0].a, 0.41518585e00, 1e-8);
        assert_delta!(bcl3.conductivities[0].b, -0.48149960e03, 1e-8);
        assert_delta!(bcl3.conductivities[0].c, 0.30788060e05, 1e-8);
        assert_delta!(bcl3.conductivities[0].d, 0.33168239e01, 1e-8);

        assert_delta!(bcl3.conductivities[1].temp_range.0, 1000.0, 1e-1);
        assert_delta!(bcl3.conductivities[1].temp_range.1, 5000.0, 1e-1);
        assert_delta!(bcl3.conductivities[1].a, 0.61148589e00, 1e-8);
        assert_delta!(bcl3.conductivities[1].b, -0.18167042e03, 1e-8);
        assert_delta!(bcl3.conductivities[1].c, -0.20976969e05, 1e-8);
        assert_delta!(bcl3.conductivities[1].d, 0.17127671e01, 1e-8);
    }
}
