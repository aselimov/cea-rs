pub fn parse_fields(line: &str, widths: &[usize]) -> Vec<String> {
    let mut fields = Vec::new();
    let mut pos = 0;

    for &width in widths {
        if let Some(field) = line.get(pos..pos + width) {
            // The replace changes the fortran formatted D exponential for the normal E exponential
            fields.push(field.trim().replace("D", "E").replace("E ", "E"));
        }
        pos += width;
    }

    fields
}
