use super::thermo_fit::{Phase, SpeciesElement, SpeciesThermoData, ThermoPolynomial};

pub fn simple_thermo_polynomial(
    cp_over_r: f64,
    h_b: f64,
    s_b: f64,
    temp_range: (f64, f64),
) -> ThermoPolynomial {
    ThermoPolynomial {
        a: vec![0.0, 0.0, cp_over_r, 0.0, 0.0, 0.0, 0.0, h_b, s_b],
        temp_range,
    }
}

pub fn h2_o2_h2o_thermo_data() -> Vec<SpeciesThermoData> {
    vec![
        SpeciesThermoData::new(
            "H2",
            vec![SpeciesElement { element: "H".to_string(), count: 2.0 }],
            Phase::Gas,
            vec![simple_thermo_polynomial(3.5, -1012.0, 29.11, (200.0, 6000.0))],
            2.01594,
            0.0,
        ),
        SpeciesThermoData::new(
            "O2",
            vec![SpeciesElement { element: "O".to_string(), count: 2.0 }],
            Phase::Gas,
            vec![simple_thermo_polynomial(3.5, -1005.0, 30.03, (200.0, 6000.0))],
            31.9988,
            0.0,
        ),
        SpeciesThermoData::new(
            "H2O",
            vec![
                SpeciesElement { element: "H".to_string(), count: 2.0 },
                SpeciesElement { element: "O".to_string(), count: 1.0 },
            ],
            Phase::Gas,
            vec![simple_thermo_polynomial(4.0, -29192.0, 23.03, (200.0, 6000.0))],
            18.01528,
            -241826.0,
        ),
    ]
}
