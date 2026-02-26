use crate::compare::types::ComparisonState;
use crate::render::sections;
use crate::state::AggregatedState;

pub fn build_report_html(state: &AggregatedState, comparison: Option<&ComparisonState>) -> String {
    let mut out = String::new();
    out.push_str("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><title>kira-organelle report</title>");
    out.push_str("<style>body{font-family:Arial,sans-serif;margin:24px;color:#111;}h1{margin:0 0 16px;}h2{margin:24px 0 8px;}h3{margin:16px 0 8px;}table{border-collapse:collapse;width:100%;margin:8px 0 16px;}th,td{border:1px solid #ccc;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#f4f4f4;}section{margin-bottom:18px;}code{background:#f4f4f4;padding:1px 4px;}</style>");
    out.push_str("</head><body><h1>kira-organelle report</h1>");

    sections::render_header(state, comparison, &mut out);
    sections::render_global_summary(state, comparison, &mut out);
    sections::render_organelle_sections(state, comparison, &mut out);
    sections::render_comparison_sections(comparison, &mut out);

    let issues = match comparison {
        Some(cmp) => &cmp.issues,
        None => &state.issues,
    };
    sections::render_issues(issues, &mut out);

    sections::render_footer(
        &state.tool.version,
        comparison.is_some(),
        comparison.map(|x| x.schema.as_str()),
        &mut out,
    );

    out.push_str("</body></html>\n");
    out
}
