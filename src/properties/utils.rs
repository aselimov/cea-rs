use std::sync::LazyLock;
use regex::Regex;

// Matches Fortran D exponent (1.23D+04, 1.23D-04, 1.23D 04) and E exponent with a space
// (1.23E 04). Captures: digit before, optional sign, digit after — drops any space.
static FORTRAN_EXP: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"([0-9])[DE]([+\-]?) ?([0-9])").unwrap()
});

pub fn parse_fields(line: &str, widths: &[usize]) -> Vec<String> {
    let mut fields = Vec::new();
    let mut pos = 0;

    for &width in widths {
        if let Some(field) = line.get(pos..pos + width) {
            fields.push(FORTRAN_EXP.replace_all(field.trim(), "${1}E${2}${3}").into_owned());
        }
        pos += width;
    }

    fields
}
