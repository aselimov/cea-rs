use std::collections::HashMap;

use crate::{
    consts::P_REF,
    properties::thermo_fit::{Phase, SpeciesThermoData},
};

pub struct MixtureComponent {
    //kg-moles component/kg_mixture
    pub(crate) n: f64,
    // Coefficients
    pub(crate) coeff: Vec<f64>,
    pub(crate) s: SpeciesThermoData,
}

pub struct GasMixture {
    pub(crate) nsum: f64,
    pub(crate) gasses: Vec<MixtureComponent>,
    pub(crate) condensed: Vec<MixtureComponent>,
    pub(crate) elements: HashMap<String, usize>,
    pub(crate) binitial: Vec<f64>,
}

impl MixtureComponent {
    // µ/RT for the current mixture component at fixed temp and pressure
    // µ = Molar Gibbs free energy = H - TS
    // Therefore µ/RT = H/RT - S/R
    pub fn chem_potential_over_rt(&self, temp: f64, pressure: f64, nsum: f64) -> f64 {
        match self.s.phase {
            Phase::Gas => {
                let p = self
                    .s
                    .polynomial_at(temp)
                    .expect("Gas doesn't have a polynomial");
                p.h_over_rt(temp) - p.s_over_r(temp)
                    + (pressure / P_REF).ln()
                    + (self.n / nsum).ln()
            }
            Phase::Condensed => todo!(),
        }
    }

    pub fn entropy_over_r(&self, temp: f64, pressure: f64, nsum: f64) -> f64 {
        match self.s.phase {
            Phase::Gas => {
                let p = self
                    .s
                    .polynomial_at(temp)
                    .expect("Gas doesn't have a polynomial");
                p.s_over_r(temp) - (self.n / nsum).ln() - (pressure / P_REF).ln()
            }
            Phase::Condensed => todo!(),
        }
    }
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

        // Now separate SpeciesThermoData in gas and condensed MixtureComponents
        let (gasses, condensed): (Vec<_>, Vec<_>) = ns
            .iter()
            .zip(species.iter())
            .map(|(n, s)| {
                let mut a = vec![0.0; elements.len()];
                s.elements.iter().for_each(|e| {
                    let i = elements.get(&e.element).unwrap();
                    a[*i] += e.count;
                });
                MixtureComponent {
                    n: *n,
                    coeff: a,
                    s: s.clone(),
                }
            })
            .partition(|c| matches!(c.s.phase, Phase::Gas));

        let binitial = get_b_current(&elements, gasses.iter().chain(condensed.iter()));

        GasMixture {
            nsum,
            gasses,
            condensed,
            elements,
            binitial,
        }
    }

    pub fn get_b_current(&self) -> Vec<f64> {
        get_b_current(
            &self.elements,
            self.gasses.iter().chain(self.condensed.iter()),
        )
    }
}

// Current kilogram-atoms of element i per kg mixture
pub fn get_b_current<'a>(
    ele_map: &HashMap<String, usize>,
    components: impl IntoIterator<Item = &'a MixtureComponent>,
) -> Vec<f64> {
    let mut b = vec![0.0; ele_map.len()];

    for c in components {
        for e in c.s.elements.iter() {
            let i = *ele_map.get(&e.element).unwrap();
            b[i] += c.n * e.count;
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

        assert_delta!(gas.gasses[0].coeff[0], 2.0, 1e-12);
        assert_delta!(gas.gasses[1].coeff[0], 0.0, 1e-12);
        assert_delta!(gas.gasses[2].coeff[0], 2.0, 1e-12);
        assert_delta!(gas.gasses[0].coeff[1], 0.0, 1e-12);
        assert_delta!(gas.gasses[1].coeff[1], 2.0, 1e-12);
        assert_delta!(gas.gasses[2].coeff[1], 1.0, 1e-12);

        let expected_b = [0.06663369409370597, 0.05830448233199272];
        assert_vec_delta!(gas.binitial, expected_b, 1e-12);
    }
}
