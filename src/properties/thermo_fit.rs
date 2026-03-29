pub enum Phase {
    Gas,
    Condensed,
}

pub struct SpeciesThermoData {
    pub name: String,
    pub elements: Vec<SpeciesElement>,
    pub phase: Phase,
    pub polynomials: Vec<ThermoPolynomial>,
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
    pub fn polynomial_at(&self, temp: f64) -> &ThermoPolynomial {
        //TODO: Not the most efficient. Can refactor to pre-compute tables
        //and do 1-d linear interpolation if needed
        //
        //TODO: I Think condensed species need to be treated differently. Verify how that works in
        //the paper.
        let i_polynomial = self
            .polynomials
            .iter()
            .rposition(|polynomial| temp > polynomial.temp_range.0)
            .unwrap_or(0);
        &self.polynomials[i_polynomial]
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
            + self.a[1] * inv_temp * inv_temp.ln()
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
        -self.a[0] * inv_temp * inv_temp - self.a[1] * inv_temp
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
        assert_delta, assert_vec_delta,
        properties::thermo_fit::{Phase, SpeciesThermoData},
    };

    #[test]
    fn test_cp_over_r() {
        // At T=1: a[0]/T^2 + a[1]/T + a[2] + a[3]*T + a[4]*T^2 + a[5]*T^3 + a[6]*T^4
        // = 1 + 2 + 3 + 4 + 5 + 6 + 7 = 28
        let result = ThermoPolynomial {
            a: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            temp_range: (0.0, 0.0),
        }
        .cp_over_r(1.0);
        assert_delta!(result, 28.0, 1e-10);

        // At T=2: 4/4 + 2/2 + 1 + 2*2 + 1*4 + 1*8 + 1*16 = 35
        let result = ThermoPolynomial {
            a: vec![4.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0],
            temp_range: (0.0, 0.0),
        }
        .cp_over_r(2.0);
        assert_delta!(result, 35.0, 1e-10);
    }

    #[test]
    fn test_h_over_rt() {
        // At T=1: ln(1/T)=0, so the a[1] term vanishes
        // -a[0] + 0 + a[2] + a[3]/2 + a[4]/3 + a[5]/4 + a[6]/5 + a[7]
        // = -1 + 0 + 3 + 12/2 + 9/3 + 8/4 + 5/5 + 8 = -1+3+6+3+2+1+8 = 22
        let result = ThermoPolynomial {
            a: vec![1.0, 2.0, 3.0, 12.0, 9.0, 8.0, 5.0, 8.0],
            temp_range: (0.0, 0.0),
        }
        .h_over_rt(1.0);
        assert_delta!(result, 22.0, 1e-10);

        // At T=2: a[1]=0 so log term vanishes; a[4..6]=0 to avoid fractions
        // -4/4 + 0 + 3 + 4*2/2 + 0 + 0 + 0 + 2/2 = -1+3+4+1 = 7
        let result = ThermoPolynomial {
            a: vec![4.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 2.0],
            temp_range: (0.0, 0.0),
        }
        .h_over_rt(2.0);
        assert_delta!(result, 7.0, 1e-10);
    }

    #[test]
    fn test_s_over_r() {
        // At T=1: ln(T)=0, so a[2] term vanishes
        // -a[0] - a[1] + 0 + a[3] + a[4]/2 + a[5]/3 + a[6]/4 + a[8]
        // = -1-2+0+4+6/2+12/3+8/4+5 = -1-2+4+3+4+2+5 = 15
        let result = ThermoPolynomial {
            a: vec![1.0, 2.0, 3.0, 4.0, 6.0, 12.0, 8.0, 0.0, 5.0],
            temp_range: (0.0, 0.0),
        }
        .s_over_r(1.0);
        assert_delta!(result, 15.0, 1e-10);

        // At T=2: a[2]=0 so log term vanishes; a[5..6]=0 to avoid fractions
        // -4/4 - 2/2 + 0 + 2*2 + 2*4/2 + 0 + 0 + 3 = -1-1+4+4+3 = 9
        let result = ThermoPolynomial {
            a: vec![4.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 3.0],
            temp_range: (0.0, 0.0),
        }
        .s_over_r(2.0);
        assert_delta!(result, 9.0, 1e-10);
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

        assert!(std::ptr::eq(data.polynomial_at(0.5), &data.polynomials[0]));
        assert!(std::ptr::eq(data.polynomial_at(50.0), &data.polynomials[0]));
        assert!(std::ptr::eq(
            data.polynomial_at(100.0),
            &data.polynomials[0]
        ));
        assert!(std::ptr::eq(
            data.polynomial_at(100.0 + 1e-12),
            &data.polynomials[1]
        ));
        assert!(std::ptr::eq(
            data.polynomial_at(150.0 + 1e-12),
            &data.polynomials[1]
        ));
        assert!(std::ptr::eq(
            data.polynomial_at(500.0 + 1e-12),
            &data.polynomials[1]
        ));
    }
}
