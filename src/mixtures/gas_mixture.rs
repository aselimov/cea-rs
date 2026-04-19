use std::collections::HashMap;

use crate::{
    consts::P_REF,
    matrix::Matrix,
    properties::thermo_fit::{Phase, SpeciesThermoData},
};

pub struct MixtureComponent {
    //kg-moles component/kg_mixture
    n: f64,
    s: SpeciesThermoData,
}

pub struct GasMixture {
    pub(crate) nsum: f64,
    pub(crate) gasses: Vec<MixtureComponent>,
    pub(crate) condensed: Vec<MixtureComponent>,
    pub(crate) coeffs: Matrix<f64>,
    pub(crate) elements: HashMap<String, usize>,
    pub(crate) binitial: Vec<f64>,
}

impl GasMixture {
    pub fn new(ns_kg_moles: &[f64], species: &[SpeciesThermoData]) -> Self {
        // First calculate the total and per species kg-moles/kg mixture
        let (nsum_kg_moles, mass_sum) = ns_kg_moles
            .iter()
            .zip(species.iter())
            .fold((0.0, 0.0), |acc, (n, s)| {
                (acc.0 + n, acc.1 + n * s.molecular_weight)
            });

        let nsum = nsum_kg_moles / mass_sum;
        let ns: Vec<f64> = ns_kg_moles.iter().map(|n| n / mass_sum).collect();

        // Now build out the element list
        let mut elements = HashMap::new();
        let mut ei = 0;
        for s in species.iter() {
            for element in s.elements.iter() {
                if !elements.contains_key(&element.element) {
                    elements.insert(element.element.clone(), ei);
                    ei += 1
                }
            }
        }

        // Now build the coefficients
        let mut coeffs = Matrix::new(ei, species.len(), 0.0);
        for (j, s) in species.iter().enumerate() {
            s.elements.iter().for_each(|e| {
                // Safe to unwrap because elements should always have e
                let i = *elements.get(&e.element).unwrap();
                coeffs.set(i, j, e.count);
            });
        }

        let binitial = get_b_current(&elements, species, &ns);

        // Now separate SpeciesThermoData in gas and condensed MixtureComponents
        let gasses = ns
            .iter()
            .zip(species.iter())
            .filter_map(|(n, s)| match s.phase {
                Phase::Gas => Some(MixtureComponent {
                    n: *n,
                    s: s.clone(),
                }),
                Phase::Condensed => None,
            })
            .collect();

        let condensed = ns
            .iter()
            .zip(species.iter())
            .filter_map(|(n, s)| match s.phase {
                Phase::Gas => None,
                Phase::Condensed => Some(MixtureComponent {
                    n: *n,
                    s: s.clone(),
                }),
            })
            .collect();

        GasMixture {
            nsum,
            gasses,
            condensed,
            coeffs,
            elements,
            binitial,
        }
    }
    // Calculate the normalized chemical potential (μ/RT) for each component in the mixture.
    // Equations 2.11 from reference paper
    pub fn gas_chem_potentials_over_rt(&self, temp: f64, pressure: f64) -> Vec<f64> {
        self.gasses
            .iter()
            .chain(self.condensed.iter())
            .map(|c| -> f64 {
                match c.s.phase {
                    Phase::Gas => {
                        let p =
                            c.s.polynomial_at(temp)
                                .expect("Gas doesn't have a polynomial");
                        p.h_over_rt(temp) - p.s_over_r(temp)
                            + (pressure / P_REF).ln()
                            + (c.n / self.nsum).ln()
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
        self.gasses
            .iter()
            .chain(self.condensed.iter())
            .map(|c| -> f64 {
                match c.s.phase {
                    Phase::Gas => {
                        let p =
                            c.s.polynomial_at(temp)
                                .expect("Gas doesn't have a polynomial");
                        p.s_over_r(temp) - (c.n / self.nsum).ln() - (pressure / P_REF).ln()
                    }
                    Phase::Condensed => todo!(),
                }
            })
            .collect()
    }

    // Calculate the normalized mixture enthalpy (H/RT)
    // Note that the enthalpy doesn't have a dependence on the pressure.
    // Equation 2.14 from the paper
    pub fn mixture_h_over_rt(&self, temp: f64) -> Vec<f64> {
        self.gasses
            .iter()
            .chain(self.condensed.iter())
            .map(|c| -> f64 {
                match c.s.phase {
                    Phase::Gas => {
                        let p =
                            c.s.polynomial_at(temp)
                                .expect("Gas doesn't have a polynomial");
                        c.n * p.h_over_rt(temp)
                    }
                    Phase::Condensed => todo!(),
                }
            })
            .collect()
    }
}

// Current kilogram-atoms of element i per kg mixture
pub fn get_b_current(
    ele_map: &HashMap<String, usize>,
    species: &[SpeciesThermoData],
    ns: &[f64],
) -> Vec<f64> {
    let mut b = vec![0.0; ele_map.len()];

    for (n, s) in ns.iter().zip(species) {
        for e in s.elements.iter() {
            let i = *ele_map.get(&e.element).unwrap();
            b[i] += n * e.count;
        }
    }
    b
}

#[cfg(test)]
mod test {
    use crate::{
        assert_delta, assert_vec_delta, mixtures::gas_mixture::GasMixture,
        properties::test_helpers::h2_o2_h2o_thermo_data,
    };

    #[test]
    pub fn test_gas_mixture_creation() {
        let species = h2_o2_h2o_thermo_data();

        let ns = [1.0, 2.0, 3.0];

        let gas = GasMixture::new(&ns, &species);

        let expected_ns = [
            0.008329211761713246,
            0.01665842352342649,
            0.02498763528513974,
        ];
        let ns: Vec<f64> = gas.gasses.iter().map(|c| c.n).collect();

        assert_vec_delta!(expected_ns, ns, 1e-12);
        assert_delta!(gas.nsum, 0.04997527057027948, 1e-12);

        assert_delta!(gas.coefs.get(0, 0).unwrap(), 2.0, 1e-12);
        assert_delta!(gas.coeffs.get(0, 1).unwrap(), 0.0, 1e-12);
        assert_delta!(gas.coeffs.get(0, 2).unwrap(), 2.0, 1e-12);
        assert_delta!(gas.coeffs.get(1, 0).unwrap(), 0.0, 1e-12);
        assert_delta!(gas.coeffs.get(1, 1).unwrap(), 2.0, 1e-12);
        assert_delta!(gas.coeffs.get(1, 2).unwrap(), 1.0, 1e-12);

        let expected_b = [0.06663369409370597, 0.05830448233199272];
        assert_vec_delta!(gas.binitial, expected_b, 1e-12);
    }
}
