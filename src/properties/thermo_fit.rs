pub enum Phase {
    Gas,
    Condensed,
}

pub struct SpeciesThermoData {
    pub name: String,
    pub elements: Vec<SpeciesElement>,
    pub phase: Phase,
    polynomials: Vec<ThermoPolynomial>,
    pub molecular_weight: f64,
    pub h_formation: f64,
}
pub struct ThermoPolynomial {
    pub a: Vec<f64>,
    pub temp_range: (f64, f64),
}

pub struct SpeciesElement {
    pub element: String,
    pub count: f64,
}

impl SpeciesThermoData {
    pub fn new(
        name: &str,
        elements: Vec<SpeciesElement>,
        phase: Phase,
        polynomials: Vec<ThermoPolynomial>,
        molecular_weight: f64,
        h_formation: f64,
    ) -> Self {
        Self {
            name: name.to_string(),
            elements,
            phase,
            polynomials,
            molecular_weight,
            h_formation,
        }
    }

    pub fn num_polynomials(&self) -> usize {
        self.polynomials.len()
    }

    pub fn polynomial_at(&self, temp: f64) -> Option<&ThermoPolynomial> {
        //TODO: Not the most efficient. Can refactor to pre-compute tables
        //and do 1-d linear interpolation if needed
        //
        //TODO: I Think condensed species need to be treated differently. Verify how that works in
        //the paper.
        if self.polynomials.is_empty() {
            return None;
        }

        let i_polynomial = self
            .polynomials
            .iter()
            .rposition(|polynomial| temp > polynomial.temp_range.0)
            .unwrap_or(0);
        Some(&self.polynomials[i_polynomial])
    }
}

impl ThermoPolynomial {
    /// Calculate using eq 4.9 from reference paper
    /// NOTE: This is normalized and unitless
    pub fn cp_over_r(&self, temp: f64) -> f64 {
        let inv_temp = 1.0 / temp;
        self.a[0] * inv_temp * inv_temp
            + self.a[1] * inv_temp
            + self.a[2]
            + self.a[3] * temp
            + self.a[4] * temp * temp
            + self.a[5] * temp * temp * temp
            + self.a[6] * temp * temp * temp * temp
    }
    /// Calculate using eq 4.10 from reference paper
    /// NOTE: This is normalized and unitless
    pub fn h_over_rt(&self, temp: f64) -> f64 {
        let inv_temp = 1.0 / temp;
        -self.a[0] * inv_temp * inv_temp
            + self.a[1] * inv_temp * temp.ln()
            + self.a[2]
            + self.a[3] * temp / 2.0
            + self.a[4] * temp * temp / 3.0
            + self.a[5] * temp * temp * temp / 4.0
            + self.a[6] * temp * temp * temp * temp / 5.0
            + self.a[7] * inv_temp
    }
    /// Calculate using eq 4.11 from reference paper
    /// NOTE: This is normalized and unitless
    pub fn s_over_r(&self, temp: f64) -> f64 {
        let inv_temp = 1.0 / temp;
        -self.a[0] * inv_temp * inv_temp / 2.0 - self.a[1] * inv_temp
            + self.a[2] * temp.ln()
            + self.a[3] * temp
            + self.a[4] * temp * temp / 2.0
            + self.a[5] * temp * temp * temp / 3.0
            + self.a[6] * temp * temp * temp * temp / 4.0
            + self.a[8]
    }
}

#[cfg(test)]
mod test {
    use super::ThermoPolynomial;
    use crate::{
        assert_delta,
        properties::thermo_fit::{Phase, SpeciesThermoData},
    };

    #[test]
    fn test_cp_over_r() {
        let result = ThermoPolynomial {
            a: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            temp_range: (0.0, 0.0),
        }
        .cp_over_r(100.0);
        assert_delta!(result, 706050403.0201, 1e-4);

        let result = ThermoPolynomial {
            a: vec![4.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0],
            temp_range: (0.0, 0.0),
        }
        .cp_over_r(2.0);
        assert_delta!(result, 35.0, 1e-10);
    }

    #[test]
    fn test_h_over_rt() {
        let result = ThermoPolynomial {
            a: vec![1.0, 2.0, 3.0, 12.0, 9.0, 8.0, 5.0, 8.0],
            temp_range: (0.0, 0.0),
        }
        .h_over_rt(2.0);
        assert_delta!(result, 63.44314718055995, 1e-10);

        let result = ThermoPolynomial {
            a: vec![4.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 2.0],
            temp_range: (0.0, 0.0),
        }
        .h_over_rt(100.0);
        assert_delta!(result, 203.0196, 1e-4);
    }

    #[test]
    fn test_s_over_r() {
        let result = ThermoPolynomial {
            a: vec![2.0, 3.0, 4.0, 1.0, 5.0, 2.0, 8.0, 1.0, 12.0],
            temp_range: (0.0, 0.0),
        }
        .s_over_r(100.0);
        assert_delta!(result, 200691797.0572474, 1e-7);

        let result = ThermoPolynomial {
            a: vec![4.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 3.0],
            temp_range: (0.0, 0.0),
        }
        .s_over_r(2.0);
        assert_delta!(result, 9.5, 1e-10);
    }

    #[test]
    fn test_polynomial_at() {
        let polynomials = vec![
            ThermoPolynomial {
                a: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
                temp_range: (1.0, 100.0),
            },
            ThermoPolynomial {
                a: vec![4.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0],
                temp_range: (100.0, 300.0),
            },
        ];

        let data = SpeciesThermoData {
            name: "".to_string(),
            elements: vec![],
            phase: Phase::Gas,
            polynomials,
            molecular_weight: 1.0,
            h_formation: 1.0,
        };

        assert!(std::ptr::eq(
            data.polynomial_at(0.5).unwrap(),
            &data.polynomials[0]
        ));
        assert!(std::ptr::eq(
            data.polynomial_at(50.0).unwrap(),
            &data.polynomials[0]
        ));
        assert!(std::ptr::eq(
            data.polynomial_at(100.0).unwrap(),
            &data.polynomials[0]
        ));
        assert!(std::ptr::eq(
            data.polynomial_at(100.0 + 1e-12).unwrap(),
            &data.polynomials[1]
        ));
        assert!(std::ptr::eq(
            data.polynomial_at(150.0 + 1e-12).unwrap(),
            &data.polynomials[1]
        ));
        assert!(std::ptr::eq(
            data.polynomial_at(500.0 + 1e-12).unwrap(),
            &data.polynomials[1]
        ));
    }
}
