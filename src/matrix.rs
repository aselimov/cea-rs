use std::fmt;
use std::{
    fmt::Display,
    ops::{Add, Div, Mul, Sub},
};

/// Numeric trait as an alias to make the generics a little bit cleaner
pub trait Numeric:
    Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Sized
    + Clone
    + Display
{
}

/// Blanket implementation for numeric types
impl<T> Numeric for T where
    T: Add<Output = Self>
        + Sub<Output = Self>
        + Mul<Output = Self>
        + Div<Output = Self>
        + Sized
        + Clone
        + Display
{
}

// Basic Error handling for Matrix operations
#[derive(Debug)]
pub enum MatrixError {
    IndexError(usize, usize, usize, usize),
}

impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MatrixError::IndexError(i, j, m, n) => write!(
                f,
                "Error accessing index [{i},{j}] for matrix with dimensions [{}, {}]",
                m, n
            ),
        }
    }
}

fn make_index_error<T>(i: usize, j: usize, m: &Matrix<T>) -> MatrixError {
    MatrixError::IndexError(i, j, m.m, m.n)
}

pub struct Matrix<T> {
    data: Vec<T>,
    m: usize,
    n: usize,
}

impl<T: Numeric> Matrix<T> {
    pub fn new(m: usize, n: usize, init_val: T) -> Self {
        let data = vec![init_val; n * m];
        Matrix { data, m, n }
    }

    fn index(&self, i: usize, j: usize) -> usize {
        i * self.n + j
    }
    pub fn get(&self, i: usize, j: usize) -> Result<&T, MatrixError> {
        self.data
            .get(self.index(i, j))
            .ok_or(make_index_error(i, j, &self))
    }

    pub fn set(&mut self, i: usize, j: usize, x: T) {
        let index = self.index(i, j);
        self.data[index] = x;
    }
}

impl<T: Numeric> Display for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let msg = (0..self.m)
            .flat_map(|i| {
                let mut parts = vec!["|".to_string()];
                parts.extend((0..self.n).map(|j| -> String {
                    if let Ok(val) = self.get(i, j) {
                        format!("{}", val)
                    } else {
                        "err".to_string()
                    }
                }));
                parts.push("|\n".to_string());
                parts
            })
            .collect::<Vec<String>>()
            .join(" ");
        write!(f, "\n{msg}")
    }
}

#[cfg(test)]
mod test {
    use crate::{assert_delta, matrix::Matrix};
    fn gen_test_matrix() -> Matrix<f64> {
        let mut m = Matrix::new(3, 4, 0.0);
        m.set(0, 0, 1.0);
        m.set(1, 1, 1.0);
        m.set(2, 2, 1.0);
        m.set(0, 3, 2.0);
        m.set(1, 3, 2.0);
        m.set(2, 3, 2.0);
        m
    }

    #[test]
    fn test_matrix_basics() {
        let m = gen_test_matrix();

        // Validate that the matrix type returns the right values
        assert_delta!(m.get(0, 0).unwrap(), 1.0, 1e-12);
        assert_delta!(m.get(1, 1).unwrap(), 1.0, 1e-12);
        assert_delta!(m.get(2, 2).unwrap(), 1.0, 1e-12);
        assert_delta!(m.get(0, 3).unwrap(), 2.0, 1e-12);
        assert_delta!(m.get(1, 3).unwrap(), 2.0, 1e-12);
        assert_delta!(m.get(2, 3).unwrap(), 2.0, 1e-12);

        // Validate that the memory is laid out as expected
        assert_delta!(m.data[0], 1.0, 1e-12);
        assert_delta!(m.data[3], 2.0, 1e-12);
        assert_delta!(m.data[5], 1.0, 1e-12);
        assert_delta!(m.data[7], 2.0, 1e-12);
        assert_delta!(m.data[10], 1.0, 1e-12);
        assert_delta!(m.data[11], 2.0, 1e-12);
    }
}
