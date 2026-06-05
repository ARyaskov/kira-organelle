use std::fmt::Write as _;

use crate::interpret::types::{
    InterpretationState, direction_str, metric_str, severity_str, source_str,
};
use crate::model::organelle::OrganelleId;
use crate::state::{AggregatedState, OrganelleState};

struct ToolGoal {
    tool: &'static str,
    goal: &'static str,
}

const CORE_TOOL_GOALS: [ToolGoal; 8] = [
    ToolGoal {
        tool: "kira-mitoqc",
        goal: "Quantify mitochondrial stress, respiratory imbalance, and decay-related quality signals.",
    },
    ToolGoal {
        tool: "kira-nuclearqc",
        goal: "Assess nucleus-associated quality stability and chromatin-linked stress metrics.",
    },
    ToolGoal {
        tool: "kira-spliceqc",
        goal: "Measure spliceosome and RNA processing quality with splicing stress indicators.",
    },
    ToolGoal {
        tool: "kira-proteoqc",
        goal: "Evaluate proteostasis load and protein quality-control pressure.",
    },
    ToolGoal {
        tool: "kira-autolys",
        goal: "Characterize autophagy/autolysosome stress and degradative compensation signals.",
    },
    ToolGoal {
        tool: "kira-secretion",
        goal: "Profile ER-Golgi and secretion pathway pressure and secretory stress state.",
    },
    ToolGoal {
        tool: "kira-energetics",
        goal: "Model energetic regimes, barrier landscape, and stability/escape structure.",
    },
    ToolGoal {
        tool: "kira-microenvironment",
        goal: "Quantify ligand-receptor microenvironment communication and group-level signaling context.",
    },
];

const OPTIONAL_TOOL_GOALS: [ToolGoal; 1] = [ToolGoal {
    tool: "kira-riboqc",
    goal: "Profile ribosome/translation-associated quality and translational burden.",
}];

pub fn build_llm_report_markdown(
    state: &AggregatedState,
    interpretation: &InterpretationState,
) -> String {
    let mut out = String::with_capacity(4096);
    out.push_str("# Kira Organelle QC Summary (LLM-Ready)\n\n");
    out.push_str("## Task\n");
    out.push_str("Aggregate organelle-oriented QC outputs into a compact narrative that an LLM can expand into a human-readable biological interpretation.\n\n");
    out.push_str("## Inputs\n");
    let _ = writeln!(&mut out, "- Primary input: `{}`", state.inputs.a);
    let _ = writeln!(
        &mut out,
        "- Comparison input: `{}`",
        state.inputs.b.as_deref().unwrap_or("-")
    );
    let _ = writeln!(&mut out, "- Timestamp: `{}`\n", state.timestamp);

    out.push_str("## Utility Goals\n");
    for item in CORE_TOOL_GOALS {
        let _ = writeln!(&mut out, "- `{}`: {}", item.tool, item.goal);
    }
    let has_ribo = state
        .organelle_states
        .iter()
        .any(|s| s.tool == "kira-riboqc");
    if has_ribo {
        for item in OPTIONAL_TOOL_GOALS {
            let _ = writeln!(&mut out, "- `{}`: {}", item.tool, item.goal);
        }
    }
    out.push('\n');

    out.push_str("## Utility Results\n");
    for item in CORE_TOOL_GOALS {
        append_tool_result(&mut out, item, state);
    }
    if has_ribo {
        for item in OPTIONAL_TOOL_GOALS {
            append_tool_result(&mut out, item, state);
        }
    }
    out.push('\n');

    out.push_str("## Cross-Utility Interpretation\n");
    let _ = writeln!(
        &mut out,
        "- Confidence score: `{:.3}`",
        interpretation.confidence
    );
    let _ = writeln!(
        &mut out,
        "- Signals detected: `{}`",
        interpretation.signals.len()
    );
    for signal in &interpretation.signals {
        let organelles = signal
            .organelles
            .iter()
            .map(|o| organelle_label(*o))
            .collect::<Vec<_>>()
            .join(", ");
        let _ = writeln!(
            &mut out,
            "- `{}`: severity={}, direction={}, organelles={}",
            signal.id,
            severity_str(signal.severity),
            direction_str(signal.direction),
            organelles
        );
        for evidence in signal.evidence.iter().take(2) {
            let _ = writeln!(
                &mut out,
                "  - evidence: source={}, organelle={}, axis={}, metric={}, value={:.3}",
                source_str(evidence.source),
                organelle_label(evidence.organelle),
                evidence.axis,
                metric_str(evidence.metric),
                evidence.value
            );
        }
    }
    if interpretation.signals.is_empty() {
        out.push_str("- No rule-based stress signals crossed alert thresholds.\n");
    }
    out.push('\n');

    out.push_str("## Notes For Downstream LLM\n");
    out.push_str(
        "- Use per-utility medians and dominant regimes to describe strongest stress drivers.\n",
    );
    out.push_str("- Mention missing utilities or contract issues as data completeness caveats.\n");
    out.push_str("- Keep conclusions descriptive and avoid clinical recommendations.\n");
    out
}

fn append_tool_result(out: &mut String, item: ToolGoal, state: &AggregatedState) {
    let _ = writeln!(out, "### {}", item.tool);
    match state.organelle_states.iter().find(|s| s.tool == item.tool) {
        Some(organelle) => {
            let _ = writeln!(
                out,
                "- Organelle: `{}`",
                organelle_label(organelle.organelle)
            );
            out.push_str("- Status: `summary loaded`\n");
            append_axes(out, organelle);
            append_regime(out, organelle);
            append_qc(out, organelle);
        }
        None => {
            out.push_str("- Status: `missing in input artifacts`\n");
            out.push_str("- Result: no summary-based metrics available for this utility.\n");
        }
    }

    let issue_count = state
        .issues
        .iter()
        .filter(|i| i.tool.as_deref() == Some(item.tool))
        .count();
    let _ = writeln!(out, "- Issues: `{issue_count}`\n");
}

fn append_axes(out: &mut String, organelle: &OrganelleState) {
    if organelle.axes.is_empty() {
        out.push_str("- Axis medians: not available.\n");
        return;
    }
    out.push_str("- Axis medians (first 3): ");
    let mut first = true;
    for axis in organelle.axes.iter().take(3) {
        if !first {
            out.push_str(", ");
        }
        first = false;
        let _ = write!(out, "{}={:.3}", axis.name, axis.median);
    }
    out.push('\n');
}

fn append_regime(out: &mut String, organelle: &OrganelleState) {
    let Some(regimes) = &organelle.regimes else {
        out.push_str("- Dominant regime: not available.\n");
        return;
    };
    let dominant = regimes.fractions.iter().max_by(|(ka, va), (kb, vb)| {
        va.partial_cmp(vb)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| ka.cmp(kb))
    });
    match dominant {
        Some((name, frac)) => {
            let _ = writeln!(out, "- Dominant regime: {name} ({frac:.3})");
        }
        None => out.push_str("- Dominant regime: not available.\n"),
    }
}

fn append_qc(out: &mut String, organelle: &OrganelleState) {
    if organelle.qc.is_empty() {
        out.push_str("- QC fractions: not available.\n");
        return;
    }
    out.push_str("- QC fractions (first 3): ");
    let mut first = true;
    for (k, v) in organelle.qc.iter().take(3) {
        if !first {
            out.push_str(", ");
        }
        first = false;
        let _ = write!(out, "{k}={v:.3}");
    }
    out.push('\n');
}

fn organelle_label(organelle: OrganelleId) -> &'static str {
    match organelle {
        OrganelleId::Mitochondria => "mitochondria",
        OrganelleId::Nucleus => "nucleus",
        OrganelleId::Spliceosome => "spliceosome",
        OrganelleId::Ribosome => "ribosome",
        OrganelleId::Proteostasis => "proteostasis",
        OrganelleId::Autophagy => "autophagy",
        OrganelleId::Secretion => "secretion",
        OrganelleId::Energetics => "energetics",
        OrganelleId::Microenvironment => "microenvironment",
    }
}
