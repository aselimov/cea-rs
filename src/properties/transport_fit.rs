pub struct SpeciesTransportData {
    pub name: String,
    pub viscosities: Vec<TransportFit>,
    pub conductivities: Vec<TransportFit>,
}

pub struct TransportFit {
    pub temp_range: (f64, f64),
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}
