use std::fmt;

use crate::{error::CEAError, matrix::Matrix};

// Basic Error handling for solver
#[derive(Debug)]
pub enum SolverError {
    NotSquareMatrix(usize, usize),
    NotColumnVector(usize, usize),
    DimensionMismatch(usize, usize),
    SingularMatrix,
}

impl fmt::Display for SolverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolverError::NotSquareMatrix(m, n) => {
                write!(f, "Matrix with dimensions [{m}, {n}] is not square")
            }
            SolverError::NotColumnVector(m, n) => write!(
                f,
                "Matrix with dimensions [{m}, {n}] is not a column vector"
            ),
            SolverError::DimensionMismatch(a_m, b_n) => write!(
                f,
                "Square matrix dimensions [{a_m}, {a_m}]\
                 do not match column vector dimensions [{b_n}, 1]"
            ),
            SolverError::SingularMatrix => write!(f, "Matrix is singular, solution does not exist"),
        }
    }
}
impl std::error::Error for SolverError {}

// Solve Ax=b for x using Gauss-Jordan elimination
// A must be a square n x n matrix and b must be a column vector, i.e. n x 1
// Algorithm taken from Numerical Methods for Engineers by Ayyub and McCuen
pub fn gauss_jordan_elimination(a: &Matrix<f64>, b: &Matrix<f64>) -> Result<Matrix<f64>, CEAError> {
    validate_ax_eq_b(a, b)?;

    let mut a = a.clone();
    let mut b = b.clone();

    for i in 0..a.shape().0 {
        // First pivot rows
        let mut max_and_index = (*a.get(i, i)?, i);
        for j in i + 1..a.shape().0 {
            max_and_index = if a.get(j, i)?.abs() > max_and_index.0.abs() {
                (*a.get(j, i)?, j)
            } else {
                max_and_index
            }
        }

        if max_and_index.0.abs() < 1e-12 {
            return Err(SolverError::SingularMatrix.into());
        }
        a.swap_rows(i, max_and_index.1);
        b.swap_rows(i, max_and_index.1);

        // Normalization
        for j in i + 1..a.shape().0 {
            a.set(i, j, a.get(i, j)? / dbg!(a.get(i, i)?));
        }
        b.set(i, 0, b.get(i, 0)? / a.get(i, i)?);
        a.set(i, i, 1.0);

        //forward pass (seems wrong)
        for k in i + 1..a.shape().0 {
            let aki = *a.get(k, i)?;
            for j in i..a.shape().0 {
                a.set(k, j, a.get(k, j)? - aki * a.get(i, j)?);
            }
            b.set(k, 0, b.get(k, 0)? - aki * b.get(i, 0)?);
        }
    }
    //backward pass
    let mut xs = b.clone();
    for i in (0..=a.shape().0 - 1).rev() {
        for j in i + 1..a.shape().0 {
            xs.set(i, 0, xs.get(i, 0)? - a.get(i, j)? * xs.get(j, 0)?);
        }
    }

    Ok(xs)
}

#[cfg(test)]
mod tests {
    use crate::{assert_delta, matrix::Matrix, solvers::equations::gauss_jordan_elimination};

    fn col_vec(data: Vec<f64>) -> Matrix<f64> {
        let n = data.len();
        Matrix::from_vec(n, 1, data).unwrap()
    }

    // 2x2: 2x + y = 5, x + 3y = 10  =>  x=1, y=3
    #[test]
    fn test_2x2_simple() {
        let a = Matrix::from_vec(2, 2, vec![2.0, 1.0, 1.0, 3.0]).unwrap();
        let b = col_vec(vec![5.0, 10.0]);
        let x = gauss_jordan_elimination(&a, &b).unwrap();
        assert_delta!(x.get(0, 0).unwrap(), 1.0, 1e-10);
        assert_delta!(x.get(1, 0).unwrap(), 3.0, 1e-10);
    }

    // 3x3 identity: Ix = b  =>  x = b
    #[test]
    fn test_3x3_identity() {
        let a = Matrix::from_vec(3, 3, vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap();
        let b = col_vec(vec![4.0, 7.0, 2.0]);
        let x = gauss_jordan_elimination(&a, &b).unwrap();
        assert_delta!(x.get(0, 0).unwrap(), 4.0, 1e-10);
        assert_delta!(x.get(1, 0).unwrap(), 7.0, 1e-10);
        assert_delta!(x.get(2, 0).unwrap(), 2.0, 1e-10);
    }

    // 3x3 general system with singularity
    // A = [[1,1,1],[0,2,1],[2,0,1]], b = [6, 7, 5]
    #[test]
    fn test_3x3_singular() {
        let a = Matrix::from_vec(3, 3, vec![1.0, 1.0, 1.0, 0.0, 2.0, 1.0, 2.0, 0.0, 1.0]).unwrap();
        let b = col_vec(vec![6.0, 7.0, 5.0]);
        assert!(gauss_jordan_elimination(&a, &b).is_err());
    }

    // 3x3 example. Example 5-6 in the reference book
    #[test]
    fn test_3x3_example_5_6() {
        let a = Matrix::from_vec(3, 3, vec![1.0, 3.0, 2.0, 2.0, 4.0, 3.0, 3.0, 4.0, 7.0]).unwrap();
        let b = col_vec(vec![15.0, 22.0, 39.0]);
        let x = gauss_jordan_elimination(&a, &b).unwrap();

        assert_delta!(x.get(0, 0).unwrap(), 1.0, 1e-12);
        assert_delta!(x.get(1, 0).unwrap(), 2.0, 1e-12);
        assert_delta!(x.get(2, 0).unwrap(), 4.0, 1e-12);
    }

    #[test]
    fn test_6x6_complex() {
        let a = Matrix::from_vec(
            6,
            6,
            vec![
                4.0, 3.0, 2.0, 1.0, 1.0, 2.0, 1.0, 5.0, 2.0, 3.0, 1.0, 1.0, 2.0, 1.0, 6.0, 2.0,
                3.0, 1.0, 1.0, 2.0, 1.0, 7.0, 2.0, 3.0, 3.0, 1.0, 2.0, 1.0, 5.0, 2.0, 2.0, 3.0,
                1.0, 2.0, 1.0, 6.0,
            ],
        )
        .unwrap();
        let b = col_vec(vec![37.0, 40.0, 51.0, 64.0, 52.0, 60.0]);
        let x = gauss_jordan_elimination(&a, &b).unwrap();
        assert_delta!(x.get(0, 0).unwrap(), 1.0, 1e-12);
        assert_delta!(x.get(1, 0).unwrap(), 2.0, 1e-12);
        assert_delta!(x.get(2, 0).unwrap(), 3.0, 1e-12);
        assert_delta!(x.get(3, 0).unwrap(), 4.0, 1e-12);
        assert_delta!(x.get(4, 0).unwrap(), 5.0, 1e-12);
        assert_delta!(x.get(5, 0).unwrap(), 6.0, 1e-12);
    }

    // Singular matrix should return an error
    #[test]
    fn test_singular_matrix() {
        let a = Matrix::from_vec(2, 2, vec![1.0, 2.0, 2.0, 4.0]).unwrap();
        let b = col_vec(vec![1.0, 2.0]);
        assert!(gauss_jordan_elimination(&a, &b).is_err());
    }

    // Non-square A should return an error
    #[test]
    fn test_non_square_matrix() {
        let a = Matrix::from_vec(2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).unwrap();
        let b = col_vec(vec![1.0, 2.0]);
        assert!(gauss_jordan_elimination(&a, &b).is_err());
    }

    // b with more than one column should return an error
    #[test]
    fn test_b_not_column_vector() {
        let a = Matrix::from_vec(2, 2, vec![1.0, 0.0, 0.0, 1.0]).unwrap();
        let b = Matrix::from_vec(2, 2, vec![1.0, 0.0, 0.0, 1.0]).unwrap();
        assert!(gauss_jordan_elimination(&a, &b).is_err());
    }

    // Dimension mismatch between A and b should return an error
    #[test]
    fn test_dimension_mismatch() {
        let a = Matrix::from_vec(3, 3, vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap();
        let b = col_vec(vec![1.0, 2.0]);
        assert!(gauss_jordan_elimination(&a, &b).is_err());
    }
}

fn validate_ax_eq_b(a: &Matrix<f64>, b: &Matrix<f64>) -> Result<(), CEAError> {
    let (a_m, a_n) = a.shape();
    if a_m != a_n {
        return Err(SolverError::NotSquareMatrix(a_m, a_n).into());
    }

    let (b_m, b_n) = b.shape();
    if b_n != 1 {
        return Err(SolverError::NotColumnVector(b_m, b_n).into());
    }

    if b_m != a_n {
        return Err(SolverError::DimensionMismatch(a_n, b_m).into());
    }
    Ok(())
}
