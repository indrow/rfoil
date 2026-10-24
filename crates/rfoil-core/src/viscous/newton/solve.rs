//! Block elimination of the Newton system.

use super::System;

impl System {
    /// Solves the system in place by block elimination.  Mass-coupling terms
    /// below the drop tolerance (`vaccel`, scaled by the arc length for the
    /// second and third equations) are skipped to save work.
    pub fn solve(&mut self, vaccel: f64, arc_length: f64) {
        let n = self.rows;
        let drop_tolerance = [vaccel, vaccel * 2.0 / arc_length, vaccel * 2.0 / arc_length];

        for row in 0..n {
            self.normalize_row(row);
            if row + 1 == n {
                break;
            }
            self.eliminate_lower_block(row);
            if row == self.te_upper_row {
                self.eliminate_wake_coupling(row);
            }
            if row + 2 >= n {
                continue;
            }
            // eliminate the dense column below the diagonal
            for target in row + 2..n {
                let coeff = *self.mass(target, row);
                for k in 0..3 {
                    if coeff[k].abs() <= drop_tolerance[k] {
                        continue;
                    }
                    for col in row + 1..n {
                        let pivot = self.mass(row, col)[2];
                        self.mass_mut(target, col)[k] -= coeff[k] * pivot;
                    }
                    for c in 0..2 {
                        self.rhs[target][k][c] -= coeff[k] * self.rhs[row][2][c];
                    }
                }
            }
        }

        // back substitution of the dense upper columns
        for row in (1..n).rev() {
            for c in 0..2 {
                let dmass = self.rhs[row][2][c];
                for above in (0..row).rev() {
                    let m = *self.mass(above, row);
                    for k in 0..3 {
                        self.rhs[above][k][c] -= m[k] * dmass;
                    }
                }
            }
        }
    }

    /// Reduces the diagonal block `[A | mass(row, row)]` of `row` to the identity.
    fn normalize_row(&mut self, row: usize) {
        let n = self.rows;
        let next = row + 1;

        // normalise equation 0, eliminate column 0 from equations 1 and 2
        let pivot = 1.0 / self.diag[row][0][0];
        self.diag[row][0][1] *= pivot;
        for col in row..n {
            self.mass_mut(row, col)[0] *= pivot;
        }
        self.scale_rhs(row, 0, pivot);
        for k in 1..3 {
            let f = self.diag[row][k][0];
            self.diag[row][k][1] -= f * self.diag[row][0][1];
            for col in row..n {
                let r0 = self.mass(row, col)[0];
                self.mass_mut(row, col)[k] -= f * r0;
            }
            self.subtract_rhs(row, k, 0, f);
        }

        // normalise equation 1, eliminate column 1 from equation 2
        let pivot = 1.0 / self.diag[row][1][1];
        for col in row..n {
            self.mass_mut(row, col)[1] *= pivot;
        }
        self.scale_rhs(row, 1, pivot);
        let f = self.diag[row][2][1];
        for col in row..n {
            let r1 = self.mass(row, col)[1];
            self.mass_mut(row, col)[2] -= f * r1;
        }
        self.subtract_rhs(row, 2, 1, f);

        // normalise equation 2, then back-eliminate column 2 and column 1
        let pivot = 1.0 / self.mass(row, row)[2];
        for col in next..n {
            self.mass_mut(row, col)[2] *= pivot;
        }
        self.scale_rhs(row, 2, pivot);
        let (f0, f1) = (self.mass(row, row)[0], self.mass(row, row)[1]);
        for col in next..n {
            let r2 = self.mass(row, col)[2];
            let m = self.mass_mut(row, col);
            m[0] -= f0 * r2;
            m[1] -= f1 * r2;
        }
        self.subtract_rhs(row, 0, 2, f0);
        self.subtract_rhs(row, 1, 2, f1);

        let f = self.diag[row][0][1];
        for col in next..n {
            let r1 = self.mass(row, col)[1];
            self.mass_mut(row, col)[0] -= f * r1;
        }
        self.subtract_rhs(row, 0, 1, f);
    }

    /// Eliminates the lower block `B` of the next row using the reduced `row`.
    fn eliminate_lower_block(&mut self, row: usize) {
        let next = row + 1;
        let b = self.lower[next];
        let below = *self.mass(next, row);
        for k in 0..3 {
            self.subtract_combination(next, k, row, [b[k][0], b[k][1], below[k]]);
        }
    }

    /// Eliminates the `Z` block coupling the first wake row to the upper TE row.
    fn eliminate_wake_coupling(&mut self, row: usize) {
        let target = self.first_wake_row;
        let z = self.wake_coupling;
        for k in 0..3 {
            self.subtract_combination(target, k, row, [z[k][0], z[k][1], 0.0]);
        }
    }

    /// `equation k of target −= Σ_j coeff[j]·(equation j of row)` over the
    /// columns right of `row`.
    fn subtract_combination(&mut self, target: usize, k: usize, row: usize, coeff: [f64; 3]) {
        for col in row + 1..self.rows {
            let src = *self.mass(row, col);
            self.mass_mut(target, col)[k] -= coeff[0] * src[0] + coeff[1] * src[1] + coeff[2] * src[2];
        }
        for c in 0..2 {
            let src = self.rhs[row];
            self.rhs[target][k][c] -= coeff[0] * src[0][c] + coeff[1] * src[1][c] + coeff[2] * src[2][c];
        }
    }

    fn scale_rhs(&mut self, row: usize, k: usize, f: f64) {
        for c in 0..2 {
            self.rhs[row][k][c] *= f;
        }
    }

    /// `rhs(row, k) −= f·rhs(row, from)` for both right-hand sides.
    fn subtract_rhs(&mut self, row: usize, k: usize, from: usize, f: f64) {
        for c in 0..2 {
            self.rhs[row][k][c] -= f * self.rhs[row][from][c];
        }
    }
}
