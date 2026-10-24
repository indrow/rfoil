//! Small dense linear-algebra kernels.

use std::ops::{Index, IndexMut};

/// Solves the tridiagonal system
///
/// ```text
///   d0 u0             = r0
///   l1 d1 u1          = r1
///      ..  ..  ..     = ..
///          ln dn      = rn
/// ```
///
/// in place: `rhs` is replaced by the solution, `diag` and `upper` are destroyed.
pub fn solve_tridiagonal(diag: &mut [f64], lower: &[f64], upper: &mut [f64], rhs: &mut [f64]) {
    let n = rhs.len();
    for k in 1..n {
        upper[k - 1] /= diag[k - 1];
        rhs[k - 1] /= diag[k - 1];
        diag[k] -= lower[k] * upper[k - 1];
        rhs[k] -= lower[k] * rhs[k - 1];
    }
    rhs[n - 1] /= diag[n - 1];
    for k in (0..n - 1).rev() {
        rhs[k] -= upper[k] * rhs[k + 1];
    }
}

/// Gaussian elimination with partial pivoting for a small `N × N` system.
/// `a` is destroyed; `rhs` is replaced by the solution.
pub fn solve_small<const N: usize>(a: &mut [[f64; N]; N], rhs: &mut [f64; N]) {
    for col in 0..N - 1 {
        let pivot_row =
            (col + 1..N).fold(col, |best, row| if a[row][col].abs() > a[best][col].abs() { row } else { best });

        // swap the pivot row into place, normalising it on the way
        let inv_pivot = 1.0 / a[pivot_row][col];
        a[pivot_row][col] = a[col][col];
        for l in col + 1..N {
            let normalised = a[pivot_row][l] * inv_pivot;
            a[pivot_row][l] = a[col][l];
            a[col][l] = normalised;
        }
        let normalised = rhs[pivot_row] * inv_pivot;
        rhs[pivot_row] = rhs[col];
        rhs[col] = normalised;

        for row in col + 1..N {
            let factor = a[row][col];
            for l in col + 1..N {
                a[row][l] -= factor * a[col][l];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    rhs[N - 1] /= a[N - 1][N - 1];
    for row in (0..N - 1).rev() {
        for k in row + 1..N {
            rhs[row] -= a[row][k] * rhs[k];
        }
    }
}

/// Dense matrix stored row-major.
#[derive(Clone, Debug)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self { rows, cols, data: vec![0.0; rows * cols] }
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.cols..(i + 1) * self.cols]
    }

    #[inline]
    pub fn row_mut(&mut self, i: usize) -> &mut [f64] {
        &mut self.data[i * self.cols..(i + 1) * self.cols]
    }
}

impl Index<(usize, usize)> for Matrix {
    type Output = f64;
    #[inline]
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        &self.data[i * self.cols + j]
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    #[inline]
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        &mut self.data[i * self.cols + j]
    }
}

/// LU factorisation with scaled partial pivoting (Crout's method).
#[derive(Clone, Debug)]
pub struct Lu {
    lu: Matrix,
    pivots: Vec<usize>,
}

impl Lu {
    /// Factors a square matrix.  A singular matrix produces non-finite
    /// entries rather than an error, as in the original solver.
    pub fn factor(mut a: Matrix) -> Self {
        let n = a.rows();
        let mut row_scale: Vec<f64> =
            (0..n).map(|i| 1.0 / a.row(i).iter().fold(0.0_f64, |m, v| m.max(v.abs()))).collect();
        let mut pivots = vec![0; n];

        for j in 0..n {
            // upper triangle of column j
            for i in 0..j {
                let mut sum = a[(i, j)];
                for k in 0..i {
                    sum -= a[(i, k)] * a[(k, j)];
                }
                a[(i, j)] = sum;
            }

            // lower triangle of column j, choosing the largest scaled pivot
            let mut best = 0.0;
            let mut pivot_row = j;
            for i in j..n {
                let mut sum = a[(i, j)];
                for k in 0..j {
                    sum -= a[(i, k)] * a[(k, j)];
                }
                a[(i, j)] = sum;
                let merit = row_scale[i] * sum.abs();
                if merit >= best {
                    pivot_row = i;
                    best = merit;
                }
            }

            if j != pivot_row {
                for k in 0..n {
                    let t = a[(pivot_row, k)];
                    a[(pivot_row, k)] = a[(j, k)];
                    a[(j, k)] = t;
                }
                row_scale[pivot_row] = row_scale[j];
            }
            pivots[j] = pivot_row;

            if j != n - 1 {
                let inv_pivot = 1.0 / a[(j, j)];
                for i in j + 1..n {
                    a[(i, j)] *= inv_pivot;
                }
            }
        }
        Self { lu: a, pivots }
    }

    /// Solves `A x = b` in place.
    pub fn solve(&self, b: &mut [f64]) {
        let n = self.pivots.len();
        let a = &self.lu;

        // forward substitution, skipping the leading zeros of b
        let mut first_nonzero: Option<usize> = None;
        for i in 0..n {
            let p = self.pivots[i];
            let mut sum = b[p];
            b[p] = b[i];
            if let Some(start) = first_nonzero {
                for j in start..i {
                    sum -= a[(i, j)] * b[j];
                }
            } else if sum != 0.0 {
                first_nonzero = Some(i);
            }
            b[i] = sum;
        }

        for i in (0..n).rev() {
            let mut sum = b[i];
            for j in i + 1..n {
                sum -= a[(i, j)] * b[j];
            }
            b[i] = sum / a[(i, i)];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lu_solves_dense_system() {
        let mut m = Matrix::zeros(3, 3);
        let rows = [[2.0, 1.0, -1.0], [-3.0, -1.0, 2.0], [-2.0, 1.0, 2.0]];
        for (i, row) in rows.iter().enumerate() {
            m.row_mut(i).copy_from_slice(row);
        }
        let lu = Lu::factor(m);
        let mut b = [8.0, -11.0, -3.0];
        lu.solve(&mut b);
        for (x, expected) in b.iter().zip([2.0, 3.0, -1.0]) {
            assert!((x - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn solve_small_pivots() {
        let mut a = [[0.0, 2.0, 1.0], [1.0, 1.0, 0.0], [3.0, 0.0, 1.0]];
        let mut rhs = [5.0, 3.0, 4.0];
        solve_small(&mut a, &mut rhs);
        for (x, expected) in rhs.iter().zip([1.0, 2.0, 1.0]) {
            assert!((x - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn solves_tridiagonal() {
        let mut diag = vec![2.0, 2.0, 2.0];
        let lower = vec![0.0, 1.0, 1.0];
        let mut upper = vec![1.0, 1.0, 0.0];
        let mut rhs = vec![3.0, 4.0, 3.0];
        solve_tridiagonal(&mut diag, &lower, &mut upper, &mut rhs);
        for x in rhs {
            assert!((x - 1.0).abs() < 1e-12);
        }
    }
}
