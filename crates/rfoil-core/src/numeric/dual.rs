//! Forward-mode dual numbers carrying `N` partial derivatives.
//!
//! The boundary-layer equations need exact Jacobians with respect to a
//! handful of primary variables.  Writing the residuals with [`Dual`]
//! arithmetic lets the chain rule be applied automatically instead of by
//! hand, which is where most of the bulk of the original Fortran lives.

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dual<const N: usize> {
    pub value: f64,
    /// Partial derivatives with respect to each independent variable.
    pub partials: [f64; N],
}

impl<const N: usize> Default for Dual<N> {
    fn default() -> Self {
        Self::constant(0.0)
    }
}

impl<const N: usize> Dual<N> {
    /// A quantity independent of all variables.
    #[inline]
    pub const fn constant(value: f64) -> Self {
        Self { value, partials: [0.0; N] }
    }

    /// The independent variable with index `slot`.
    #[inline]
    pub fn var(value: f64, slot: usize) -> Self {
        let mut partials = [0.0; N];
        partials[slot] = 1.0;
        Self { value, partials }
    }

    /// `f(self)` given the value `f` and the derivative `df/dself`.
    #[inline]
    fn apply(self, f: f64, df: f64) -> Self {
        Self { value: f, partials: self.partials.map(|p| p * df) }
    }

    /// Builds `f(a, b, ...)` from its value and the pairs `(∂f/∂a, a), (∂f/∂b, b), ...`.
    #[inline]
    pub fn chain(value: f64, terms: &[(f64, Self)]) -> Self {
        let mut partials = [0.0; N];
        for (df, arg) in terms {
            for (out, p) in partials.iter_mut().zip(arg.partials) {
                *out += df * p;
            }
        }
        Self { value, partials }
    }

    #[inline]
    pub fn ln(self) -> Self {
        self.apply(self.value.ln(), 1.0 / self.value)
    }

    #[inline]
    pub fn exp(self) -> Self {
        let e = self.value.exp();
        self.apply(e, e)
    }

    #[inline]
    pub fn sqrt(self) -> Self {
        let r = self.value.sqrt();
        self.apply(r, 0.5 / r)
    }

    #[inline]
    pub fn powf(self, p: f64) -> Self {
        self.apply(self.value.powf(p), p * self.value.powf(p - 1.0))
    }

    #[inline]
    pub fn tanh(self) -> Self {
        let t = self.value.tanh();
        self.apply(t, 1.0 - t * t)
    }

    /// Weighted blend `(1 − w)·self + w·other` with a dual weight.
    #[inline]
    pub fn lerp(self, other: Self, w: Self) -> Self {
        self * (1.0 - w) + other * w
    }
}

impl<const N: usize> From<f64> for Dual<N> {
    fn from(value: f64) -> Self {
        Self::constant(value)
    }
}

impl<const N: usize> Neg for Dual<N> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self { value: -self.value, partials: self.partials.map(|p| -p) }
    }
}

impl<const N: usize> Add for Dual<N> {
    type Output = Self;
    #[inline]
    fn add(mut self, rhs: Self) -> Self {
        self.value += rhs.value;
        for (a, b) in self.partials.iter_mut().zip(rhs.partials) {
            *a += b;
        }
        self
    }
}

impl<const N: usize> AddAssign for Dual<N> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const N: usize> Sub for Dual<N> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl<const N: usize> Mul for Dual<N> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let mut partials = [0.0; N];
        for i in 0..N {
            partials[i] = self.partials[i] * rhs.value + self.value * rhs.partials[i];
        }
        Self { value: self.value * rhs.value, partials }
    }
}

impl<const N: usize> Div for Dual<N> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        let inv = 1.0 / rhs.value;
        let quotient = self.value * inv;
        let mut partials = [0.0; N];
        for i in 0..N {
            partials[i] = (self.partials[i] - quotient * rhs.partials[i]) * inv;
        }
        Self { value: quotient, partials }
    }
}

/// Mixed dual/scalar arithmetic in both operand orders.
macro_rules! scalar_ops {
    ($($trait:ident $method:ident),*) => {$(
        impl<const N: usize> $trait<f64> for Dual<N> {
            type Output = Self;
            #[inline]
            fn $method(self, rhs: f64) -> Self { $trait::$method(self, Self::constant(rhs)) }
        }
        impl<const N: usize> $trait<Dual<N>> for f64 {
            type Output = Dual<N>;
            #[inline]
            fn $method(self, rhs: Dual<N>) -> Dual<N> { $trait::$method(Dual::constant(self), rhs) }
        }
    )*};
}
scalar_ops!(Add add, Sub sub, Mul mul, Div div);

#[cfg(test)]
mod tests {
    use super::*;

    type D2 = Dual<2>;

    #[test]
    fn product_and_quotient_rules() {
        let x = D2::var(3.0, 0);
        let y = D2::var(2.0, 1);
        let f = x * x / y + (x * y).ln();
        assert!((f.value - (4.5 + 6.0_f64.ln())).abs() < 1e-12);
        assert!((f.partials[0] - (3.0 + 1.0 / 3.0)).abs() < 1e-12);
        assert!((f.partials[1] - (-9.0 / 4.0 + 0.5)).abs() < 1e-12);
    }

    #[test]
    fn chain_matches_direct_composition() {
        let x = D2::var(0.7, 0);
        let t = 0.7_f64.tanh();
        assert_eq!(x.tanh(), D2::chain(t, &[(1.0 - t * t, x)]));
    }
}
