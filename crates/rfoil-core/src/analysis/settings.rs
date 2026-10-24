//! Operating-condition settings.

use crate::viscous::{Amplification, Calibration};

/// How Mach and Reynolds numbers vary with CL.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Scaling {
    /// Type 1: fixed Mach and Reynolds number.
    #[default]
    Fixed,
    /// Type 2: fixed `M·√CL` and `Re·√CL` (fixed wing loading, varying speed).
    FixedLift,
    /// Type 3: fixed Mach and `Re·CL` (fixed loading and speed, varying chord).
    FixedLiftAndSpeed,
}

/// Analysis settings.
#[derive(Clone, Debug, PartialEq)]
pub struct Settings {
    /// Chord Reynolds number; `0` selects an inviscid analysis.
    pub reynolds: f64,
    pub mach: f64,
    pub scaling: Scaling,
    /// Critical amplification exponent of the eᴺ method.
    pub ncrit: f64,
    /// Forced transition x/c on the top and bottom surfaces (`1` = free).
    pub xtr_top: f64,
    pub xtr_bottom: f64,
    /// Newton iteration limit per operating point.
    pub max_iterations: usize,
    /// Drop tolerance of the block solver.
    pub vaccel: f64,
    /// Wake length in chords.
    pub wake_length: f64,
    /// Moment reference point.
    pub moment_ref: (f64, f64),
    pub amplification: Amplification,
    pub calibration: Calibration,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            reynolds: 0.0,
            mach: 0.0,
            scaling: Scaling::Fixed,
            ncrit: 9.0,
            xtr_top: 1.0,
            xtr_bottom: 1.0,
            max_iterations: 50,
            vaccel: 0.01,
            wake_length: 1.0,
            moment_ref: (0.25, 0.0),
            amplification: Amplification::Envelope,
            calibration: Calibration::default(),
        }
    }
}

/// Mach and Reynolds numbers at a given CL, with their CL-derivatives.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowNumbers {
    pub mach: f64,
    pub mach_cl: f64,
    pub reynolds: f64,
    pub reynolds_cl: f64,
}

impl Settings {
    pub fn is_viscous(&self) -> bool {
        self.reynolds > 0.0
    }

    /// Actual Mach and Reynolds numbers at lift coefficient `cl`.
    pub fn flow_numbers(&self, cl: f64) -> FlowNumbers {
        let cl = cl.max(1.0e-6); // the scalings divide by CL
        let (mut mach, mut mach_cl) = match self.scaling {
            Scaling::FixedLift => {
                let m = self.mach / cl.sqrt();
                (m, -0.5 * m / cl)
            }
            _ => (self.mach, 0.0),
        };
        let (mut reynolds, mut reynolds_cl) = match self.scaling {
            Scaling::Fixed => (self.reynolds, 0.0),
            Scaling::FixedLift => {
                let r = self.reynolds / cl.sqrt();
                (r, -0.5 * r / cl)
            }
            Scaling::FixedLiftAndSpeed => {
                let r = self.reynolds / cl;
                (r, -r / cl)
            }
        };
        // keep the flow subsonic and the Reynolds number bounded at tiny CL
        if mach >= 0.99 {
            (mach, mach_cl) = (0.99, 0.0);
        }
        if self.reynolds > 0.0 && reynolds / self.reynolds > 100.0 {
            (reynolds, reynolds_cl) = (100.0 * self.reynolds, 0.0);
        }
        FlowNumbers { mach, mach_cl, reynolds, reynolds_cl }
    }

    /// True if changing from `self` to `other` invalidates a converged viscous solution.
    pub(crate) fn changes_viscous_problem(&self, other: &Settings) -> bool {
        self.reynolds != other.reynolds
            || self.mach != other.mach
            || self.scaling != other.scaling
            || self.ncrit != other.ncrit
            || self.xtr_top != other.xtr_top
            || self.xtr_bottom != other.xtr_bottom
            || self.amplification != other.amplification
            || self.calibration != other.calibration
            || self.moment_ref != other.moment_ref
    }
}
