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
