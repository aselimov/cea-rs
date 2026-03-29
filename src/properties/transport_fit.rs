pub struct SpeciesTransportData {
    pub name: String,
    pub viscosities: Vec<TransportFit>,
    pub conductivities: Vec<TransportFit>,
}

impl SpeciesTransportData {
    pub fn viscosity_at(&self, temp: f64) -> f64 {
        let i_viscosity = self
            .viscosities
            .iter()
            .rposition(|viscosity| temp > viscosity.temp_range.0)
            .unwrap_or(0);
        self.viscosities[i_viscosity].compute(temp)
    }

    pub fn conductivity_at(&self, temp: f64) -> f64 {
        let i_conductivity = self
            .conductivities
            .iter()
            .rposition(|conductivity| temp > conductivity.temp_range.0)
            .unwrap_or(0);
        self.conductivities[i_conductivity].compute(temp)
    }
}

pub struct TransportFit {
    pub temp_range: (f64, f64),
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}

impl TransportFit {
    pub fn compute(&self, temp: f64) -> f64 {
        let inv_temp = 1.0 / temp;
        self.a * temp.ln() + self.b * inv_temp + self.c * inv_temp * inv_temp + self.d
    }
}

#[cfg(test)]
mod test {
    use crate::{
        assert_delta,
        properties::transport_fit::{SpeciesTransportData, TransportFit},
    };

    #[test]
    fn test_transport_fit_compute() {
        let fit = TransportFit {
            temp_range: (1.0, 2.0),
            a: 10.0,
            b: 20.0,
            c: 30.0,
            d: 40.0,
        };

        assert_delta!(fit.compute(1.0), 90.0, 1e-12);
    }

    #[test]
    fn test_calc_transport_properties() {
        let viscosities = vec![
            TransportFit {
                temp_range: (1.0, 100.0),
                a: 10.0,
                b: 20.0,
                c: 30.0,
                d: 40.0,
            },
            TransportFit {
                temp_range: (100.0, 200.0),
                a: 1.0,
                b: 2.0,
                c: 3.0,
                d: 4.0,
            },
        ];

        let conductivities = vec![
            TransportFit {
                temp_range: (1.0, 100.0),
                a: 10.0,
                b: 20.0,
                c: 30.0,
                d: 40.0,
            },
            TransportFit {
                temp_range: (100.0, 200.0),
                a: 1.0,
                b: 2.0,
                c: 3.0,
                d: 4.0,
            },
        ];

        let data = SpeciesTransportData {
            viscosities,
            conductivities,
            name: "".to_string(),
        };

        assert_delta!(
            data.conductivity_at(0.5),
            data.conductivities[0].compute(0.5),
            1e-12
        );
        assert_delta!(
            data.conductivity_at(1.0),
            data.conductivities[0].compute(1.0),
            1e-12
        );
        assert_delta!(
            data.conductivity_at(50.0),
            data.conductivities[0].compute(50.0),
            1e-12
        );
        assert_delta!(
            data.conductivity_at(100.0),
            data.conductivities[0].compute(100.0),
            1e-12
        );
        assert_delta!(
            data.conductivity_at(100.0 + 1e-12),
            data.conductivities[1].compute(100.0 + 1e-12),
            1e-12
        );

        assert_delta!(
            data.conductivity_at(200.0),
            data.conductivities[1].compute(200.0),
            1e-12
        );
        assert_delta!(
            data.conductivity_at(500.0),
            data.conductivities[1].compute(500.0),
            1e-12
        );

        assert_delta!(
            data.viscosity_at(0.5),
            data.viscosities[0].compute(0.5),
            1e-12
        );
        assert_delta!(
            data.viscosity_at(1.0),
            data.viscosities[0].compute(1.0),
            1e-12
        );
        assert_delta!(
            data.viscosity_at(50.0),
            data.viscosities[0].compute(50.0),
            1e-12
        );
        assert_delta!(
            data.viscosity_at(100.0),
            data.viscosities[0].compute(100.0),
            1e-12
        );
        assert_delta!(
            data.viscosity_at(100.0 + 1e-12),
            data.viscosities[1].compute(100.0 + 1e-12),
            1e-12
        );

        assert_delta!(
            data.viscosity_at(200.0),
            data.viscosities[1].compute(200.0),
            1e-12
        );
        assert_delta!(
            data.viscosity_at(500.0),
            data.viscosities[1].compute(500.0),
            1e-12
        );
    }
}
