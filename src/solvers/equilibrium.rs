use crate::{matrix::Matrix, mixtures::gas_mixture::GasMixture};

pub enum Problem {
    TP { temp: f64, pressure: f64 },
}

const ITERATIONS: usize = 1000;
const TOLERANCE: f64 = 1e-8;

pub fn solve(gas: &GasMixture, problem: Problem) -> GasMixture {
    todo!()
}

// Construct the matrices representing the iteration equations that can be solved for variable
// updates.
//
// NOTE: This is not the most computationally efficient implementation, but it is more readable.
// Future versions may combine loops if performance is a big issue.
pub fn construct_matrices(
    gas: &GasMixture,
    problem: Problem,
    temp: f64,
    pressure: f64,
) -> (Matrix<f64>, Matrix<f64>) {
    let equations = gas.elements.len() + gas.condensed.len();

    match problem {
        Problem::TP {
            temp: _,
            pressure: _,
        } => (),
    }

    let mut a = Matrix::new(equations, equations, 0.0);
    let mut y = Matrix::new(equations, 1, 0.0);

    let gas_potentials: Vec<_> = gas
        .gasses
        .iter()
        .map(|g| g.chem_potential_over_rt(temp, pressure, gas.nsum))
        .collect();

    let condensed_potentials: Vec<_> = gas
        .gasses
        .iter()
        .map(|g| g.chem_potential_over_rt(temp, pressure, gas.nsum))
        .collect();

    let ngas = gas
        .gasses
        .iter()
        .chain(gas.condensed.iter())
        .fold(0.0, |acc, g| acc + g.n);

    // Equations 2.24
    for k in 0..gas.elements.len() {
        for i in 0..gas.elements.len() {
            for (j, g) in gas.gasses.iter().enumerate() {
                a.set(k, i, a.get(k, i).unwrap() + g.coeff[k] * g.coeff[i] * g.n);
            }
        }

        let del_nj_start = gas.elements.len();
        for (j, c) in gas.condensed.iter().enumerate() {
            a.set(k, j + del_nj_start, a.get(k, j).unwrap() + c.coeff[k])
        }

        let del_ln_n_start = del_nj_start + gas.condensed.len();
        for (j, g) in gas.gasses.iter().enumerate() {
            a.set(
                k,
                del_ln_n_start,
                a.get(k, del_ln_n_start).unwrap() + g.n * g.coeff[j],
            );
        }
        //TODO: \delta T term
        let y_k = gas
            .gasses
            .iter()
            .zip(gas_potentials.iter())
            .fold(gas.binitial[k] - gas.get_b_current()[k], |acc, (c, mu)| {
                acc + c.n * c.coeff[k] * mu
            });
        y.set(k, 0, y_k);
    }

    // Equation 2.25
    for (j, (c, mu)) in gas
        .condensed
        .iter()
        .zip(condensed_potentials.iter())
        .enumerate()
    {
        for i in 0..gas.elements.len() {
            a.set(i, j, a.get(i, j).unwrap() + c.coeff[i]);
        }
        //TODO: \delta T term
        y.set(j + gas.elements.len(), 1, *mu);
    }

    // Equation 2.26
    let offset = gas.elements.len() + gas.condensed.len();
    for (j, (g, mu)) in gas.gasses.iter().zip(gas_potentials.iter()).enumerate() {
        for i in 0..gas.elements.len() {
            a.set(offset, i, a.get(offset, i).unwrap() + g.n * g.coeff[i]);
        }
        y.set(offset, 0, y.get(offset, 0).unwrap() + g.n * mu);
    }
    let del_ln_n_start = gas.elements.len() + gas.condensed.len();
    a.set(offset, del_ln_n_start, ngas - gas.nsum);
    y.set(offset, 0, y.get(offset, 0).unwrap() + gas.nsum - ngas);

    (a, y)
}
