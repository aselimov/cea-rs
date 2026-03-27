use std::fmt;

#[derive(Debug)]
pub enum PropertiesError {
    ParseError(String, String, String),
    InvalidFormat(String),
    InvalidFile,
    InvalidLine(String),
}

impl fmt::Display for PropertiesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PropertiesError::ParseError(var, var_type, line) => {
                write!(
                    f,
                    "Failed to parse {} with type {}\nline is\n{}",
                    var, var_type, line
                )
            }
            PropertiesError::InvalidFormat(msg) => write!(f, "invalid format: {}", msg),
            PropertiesError::InvalidLine(msg) => {
                write!(f, "Failed to split line when parsing {}", msg)
            }
            PropertiesError::InvalidFile => write!(f, "Not enough lines in file"),
        }
    }
}

impl std::error::Error for PropertiesError {}

pub fn make_parse_error(var_name: &str, var_type: &str, line: &str) -> PropertiesError {
    PropertiesError::ParseError(var_name.to_string(), var_type.to_string(), line.to_string())
}
