use thiserror::Error;

/// Errors reported by the solver library.
#[derive(Debug, Error)]
pub enum Error {
    #[error("invalid NACA designation `{0}` (use 4 digits, or 5 digits starting with 210..250)")]
    InvalidNaca(String),

    #[error("airfoil file contains too few coordinate points ({0})")]
    TooFewPoints(usize),

    #[error("could not parse airfoil coordinates: {0}")]
    Parse(String),

    #[error("no airfoil has been loaded")]
    NoAirfoil,

    #[error("invalid parameter: {0}")]
    InvalidParameter(String),
}

pub type Result<T> = std::result::Result<T, Error>;
