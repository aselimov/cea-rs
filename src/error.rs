use std::fmt;

use crate::{matrix::MatrixError, properties::PropertiesError, solvers::SolverError};

#[derive(Debug)]
pub enum CEAError {
    Matrix(MatrixError),
    Solver(SolverError),
    Properties(PropertiesError),
}

impl fmt::Display for CEAError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CEAError::Matrix(e) => write!(f, "matrix error: {}", e),
            CEAError::Solver(e) => write!(f, "solver error: {}", e),
            CEAError::Properties(e) => write!(f, "properties error: {}", e),
        }
    }
}

impl std::error::Error for CEAError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            CEAError::Matrix(e) => Some(e),
            CEAError::Solver(e) => Some(e),
            CEAError::Properties(e) => Some(e),
        }
    }
}

impl From<MatrixError> for CEAError {
    fn from(e: MatrixError) -> Self {
        CEAError::Matrix(e)
    }
}

impl From<SolverError> for CEAError {
    fn from(e: SolverError) -> Self {
        CEAError::Solver(e)
    }
}

impl From<PropertiesError> for CEAError {
    fn from(e: PropertiesError) -> Self {
        CEAError::Properties(e)
    }
}
