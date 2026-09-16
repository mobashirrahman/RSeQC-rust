use thiserror::Error;

#[derive(Debug, Error)]
pub enum CoreError {
    #[error("invalid coordinate: {0}")]
    InvalidCoordinate(String),
}
