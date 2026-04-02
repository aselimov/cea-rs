use crate::{
    consts::P_REF,
    properties::{
        thermo_fit::{Phase, SpeciesThermoData},
        transport_fit::SpeciesTransportData,
    },
};

pub struct GasMixture {
    pub(crate) ns: Vec<f64>,
    pub(crate) nsum: f64
    pub(crate) species: Vec<SpeciesThermoData>,
    pub(crate) transport_data: Vec<SpeciesTransportData>,
}

impl GasMixture {
    // Calculate the normalized chemical potential (μ/RT) for each component in the mixture.
    //
    // Equations 2.11 from reference paper
    pub fn gas_chem_potentials_over_rt(&self, temp: f64, pressure: f64) -> Vec<f64> {
        self.ns
            .iter()
            .zip(self.species.iter())
            .map(|(n, s)| -> f64 {
                match s.phase {
                    Phase::Gas => {
                        let p = s
                            .polynomial_at(temp)
                            .expect("Gas doesn't have a polynomial");
                        p.h_over_rt(temp) - p.s_over_r(temp) + (pressure / P_REF).ln() + (n/self.nsum).ln()
                    }
                    Phase::Condensed => todo!(),
                }
            })
            .collect()
    }
    
    // Calculate the normalized entropy (S/R) for each mixture component
    //
    // Equations 2.17 from reference paper
    pub fn gas_entropies_over_rt(&self, temp: f64, pressure: f64) -> Vec<f64> {
        self.ns
            .iter()
            .zip(self.species.iter())
            .map(|(n, s)| -> f64 {
                match s.phase {
                    Phase::Gas => {
                        let p = s
                            .polynomial_at(temp)
                            .expect("Gas doesn't have a polynomial");
                        p.s_over_r(temp) - (n/self.nsum).ln() - (pressure/P_REF).ln()
                    }
                    Phase::Condensed => todo!(),
                }
            })
            .collect()
    }
    
    // Calculate the normalized mixture enthalpy (H/RT)
    // Note that the enthalpy doesn't have a dependence on the pressure.
    // Equation 2.14 from the paper
    pub fn mixture_h_over_rt(&self, temp: f64) -> Vec<f64>{
         self.ns
            .iter()
            .zip(self.species.iter())
            .map(|(n, s)| -> f64 {
                match s.phase {
                    Phase::Gas => {
                        let p = s
                            .polynomial_at(temp)
                            .expect("Gas doesn't have a polynomial");
                        n*p.h_over_rt(temp)
                    }
                    Phase::Condensed => todo!(),
                }
            })
            .collect()       
    }
}
