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
    + Copy
    + Display
    + Default
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
        + Default
        + Copy
{
}

// Basic Error handling for Matrix operations
#[derive(Debug)]
pub enum MatrixError {
    IndexError(usize, usize, usize, usize),
    AddError(usize, usize, usize, usize),
    MultiplicationError(usize, usize, usize, usize),
    FromVecError(usize, usize, usize),
}

impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MatrixError::IndexError(i, j, m, n) => write!(
                f,
                "Error accessing index [{i},{j}] for matrix with dimensions [{}, {}]",
                m, n
            ),

            MatrixError::MultiplicationError(i, j, m, n) => write!(
                f,
                "Matrices with dimensions [{i},{j}] and [{m},{n}] cannot be multiplied",
            ),
            MatrixError::AddError(i, j, m, n) => write!(
                f,
                "Matrices with dimensions [{i},{j}] and [{m},{n}] cannot be added",
            ),
            MatrixError::FromVecError(i, j, len) => write!(
                f,
                "Matrices with dimensions [{i},{j}] cannot be created from vec with len={len}",
            ),
        }
    }
}

impl std::error::Error for MatrixError {}

fn make_index_error<T>(i: usize, j: usize, m: &Matrix<T>) -> MatrixError {
    MatrixError::IndexError(i, j, m.m, m.n)
}

#[derive(Clone)]
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

    pub fn from_vec(m: usize, n: usize, data: Vec<T>) -> Result<Self, MatrixError> {
        if m * n != data.len() {
            return Err(MatrixError::FromVecError(m, n, data.len()));
        }
        Ok(Self { m, n, data })
    }

    fn index(&self, i: usize, j: usize) -> usize {
        i * self.n + j
    }
    pub fn get(&self, i: usize, j: usize) -> Result<&T, MatrixError> {
        if i >= self.m || j >= self.n {
            return Err(make_index_error(i, j, self));
        }
        self.data
            .get(self.index(i, j))
            .ok_or(make_index_error(i, j, self))
    }

    pub fn set(&mut self, i: usize, j: usize, x: T) {
        let index = self.index(i, j);
        self.data[index] = x;
    }

    pub fn add(&self, other: &Matrix<T>) -> Result<Matrix<T>, MatrixError> {
        // Compatibility check
        if self.m != other.m || self.n != other.n {
            return Err(MatrixError::AddError(self.m, self.n, other.m, other.n));
        }

        let data = self
            .data
            .iter()
            .zip(other.data.iter())
            .map(|(a, b)| *a + *b)
            .collect();
        Self::from_vec(self.m, self.n, data)
    }

    pub fn mul(&self, other: &Matrix<T>) -> Result<Matrix<T>, MatrixError> {
        let mut c = Matrix::new(self.m, other.n, T::default());
        // Compatibility check
        if self.n != other.m {
            return Err(MatrixError::MultiplicationError(
                self.m, self.n, other.m, other.n,
            ));
        }

        for i in 0..self.m {
            for k in 0..self.n {
                for j in 0..other.n {
                    c.set(
                        i,
                        j,
                        *c.get(i, j)? + (*self.get(i, k)?) * (*other.get(k, j)?),
                    );
                }
            }
        }
        Ok(c)
    }

    pub fn transpose(&self) -> Result<Self, MatrixError> {
        let mut c = Self::new(self.n, self.m, T::default());

        for i in 0..self.m {
            for j in 0..self.n {
                c.set(j, i, *self.get(i, j)?);
            }
        }
        Ok(c)
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.m, self.n)
    }

    pub fn swap_rows(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }
        for col in 0..self.n {
            let first_index = self.index(i, col);
            let second_index = self.index(j, col);
            self.data.swap(first_index, second_index);
        }
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
        write!(f, "\n {msg}")
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

        // Test from_vec
        let data = vec![1, 2, 3, 4, 5, 6];
        let m = Matrix::from_vec(3, 2, data).unwrap();
        assert_eq!(*m.get(0, 0).unwrap(), 1);
        assert_eq!(*m.get(0, 1).unwrap(), 2);
        assert_eq!(*m.get(1, 0).unwrap(), 3);
        assert_eq!(*m.get(1, 1).unwrap(), 4);
        assert_eq!(*m.get(2, 0).unwrap(), 5);
        assert_eq!(*m.get(2, 1).unwrap(), 6);
    }

    #[test]
    fn test_matrix_math() {
        let a = Matrix::from_vec(3, 2, vec![1, 2, 3, 4, 5, 6]).unwrap();
        let b = Matrix::from_vec(3, 2, vec![0, 1, 0, 1, 0, 1]).unwrap();

        let c = a.add(&b).unwrap();
        assert_eq!(*c.get(0, 0).unwrap(), 1);
        assert_eq!(*c.get(0, 1).unwrap(), 3);
        assert_eq!(*c.get(1, 0).unwrap(), 3);
        assert_eq!(*c.get(1, 1).unwrap(), 5);
        assert_eq!(*c.get(2, 0).unwrap(), 5);
        assert_eq!(*c.get(2, 1).unwrap(), 7);

        assert!(a.mul(&b).is_err());

        let a = Matrix::from_vec(3, 2, vec![1, 2, 3, 4, 5, 6]).unwrap();
        let b = Matrix::from_vec(2, 3, vec![0, 1, 0, 1, 0, 1]).unwrap();
        let c = a.mul(&b).unwrap();

        assert_eq!(c.data, vec![2, 1, 2, 4, 3, 4, 6, 5, 6]);

        let a = Matrix::from_vec(3, 4, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]).unwrap();
        let atranspose = a.transpose().unwrap();

        assert_eq!(*atranspose.get(0, 0).unwrap(), 1);
        assert_eq!(*atranspose.get(0, 1).unwrap(), 5);
        assert_eq!(*atranspose.get(0, 2).unwrap(), 9);
        assert!(atranspose.get(0, 3).is_err());
        assert_eq!(*atranspose.get(1, 0).unwrap(), 2);
        assert_eq!(*atranspose.get(1, 1).unwrap(), 6);
        assert_eq!(*atranspose.get(1, 2).unwrap(), 10);
        assert_eq!(*atranspose.get(2, 0).unwrap(), 3);
        assert_eq!(*atranspose.get(2, 1).unwrap(), 7);
        assert_eq!(*atranspose.get(2, 2).unwrap(), 11);
        assert_eq!(*atranspose.get(3, 0).unwrap(), 4);
        assert_eq!(*atranspose.get(3, 1).unwrap(), 8);
        assert_eq!(*atranspose.get(3, 2).unwrap(), 12);
    }

    #[test]
    fn test_swap_rows() {
        let mut m = Matrix::from_vec(3, 3, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]).unwrap();

        // Swap rows 0 and 2
        m.swap_rows(0, 2);
        assert_eq!(*m.get(0, 0).unwrap(), 7);
        assert_eq!(*m.get(0, 1).unwrap(), 8);
        assert_eq!(*m.get(0, 2).unwrap(), 9);
        assert_eq!(*m.get(2, 0).unwrap(), 1);
        assert_eq!(*m.get(2, 1).unwrap(), 2);
        assert_eq!(*m.get(2, 2).unwrap(), 3);
        // Row 1 unchanged
        assert_eq!(*m.get(1, 0).unwrap(), 4);
        assert_eq!(*m.get(1, 1).unwrap(), 5);
        assert_eq!(*m.get(1, 2).unwrap(), 6);

        // Swapping a row with itself is a no-op
        m.swap_rows(1, 1);
        assert_eq!(*m.get(1, 0).unwrap(), 4);
        assert_eq!(*m.get(1, 1).unwrap(), 5);
        assert_eq!(*m.get(1, 2).unwrap(), 6);
    }
}
