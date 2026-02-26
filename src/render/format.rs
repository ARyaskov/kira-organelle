pub fn escape_html(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for ch in input.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(ch),
        }
    }
    out
}

pub fn fmt_f64(value: f64) -> String {
    format!("{value:.3}")
}

pub fn fmt_opt_f64(value: Option<f64>) -> String {
    match value {
        Some(v) => fmt_f64(v),
        None => "-".to_string(),
    }
}

pub fn fmt_delta(value: f64) -> String {
    if value < 0.0 {
        format!("&minus;{}", fmt_f64(value.abs()))
    } else {
        fmt_f64(value)
    }
}

pub fn organelle_label(raw: &str) -> &'static str {
    match raw {
        "mitochondria" => "Mitochondria",
        "nucleus" => "Nucleus",
        "spliceosome" => "Spliceosome",
        "ribosome" => "Ribosome",
        "proteostasis" => "Proteostasis",
        "autophagy" => "Autophagy",
        "secretion" => "Secretion",
        _ => "Unknown",
    }
}
