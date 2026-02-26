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
    let mut out = String::new();
    out.push_str("# Kira Organelle QC Summary (LLM-Ready)\n\n");
    out.push_str("## Task\n");
    out.push_str("Aggregate organelle-oriented QC outputs into a compact narrative that an LLM can expand into a human-readable biological interpretation.\n\n");
    out.push_str("## Inputs\n");
    out.push_str(&format!("- Primary input: `{}`\n", state.inputs.a));
    out.push_str(&format!(
        "- Comparison input: `{}`\n",
        state.inputs.b.as_deref().unwrap_or("-")
    ));
    out.push_str(&format!("- Timestamp: `{}`\n\n", state.timestamp));

    out.push_str("## Utility Goals\n");
    for item in CORE_TOOL_GOALS {
        out.push_str(&format!("- `{}`: {}\n", item.tool, item.goal));
    }
    if state
        .organelle_states
        .iter()
        .any(|s| s.tool == "kira-riboqc")
    {
        for item in OPTIONAL_TOOL_GOALS {
            out.push_str(&format!("- `{}`: {}\n", item.tool, item.goal));
        }
    }
    out.push('\n');

    out.push_str("## Utility Results\n");
    for item in CORE_TOOL_GOALS {
        append_tool_result(&mut out, item, state);
    }
    if state
        .organelle_states
        .iter()
        .any(|s| s.tool == "kira-riboqc")
    {
        for item in OPTIONAL_TOOL_GOALS {
            append_tool_result(&mut out, item, state);
        }
    }
    out.push('\n');

    out.push_str("## Cross-Utility Interpretation\n");
    out.push_str(&format!(
        "- Confidence score: `{:.3}`\n",
        interpretation.confidence
    ));
    out.push_str(&format!(
        "- Signals detected: `{}`\n",
        interpretation.signals.len()
    ));
    for signal in &interpretation.signals {
        out.push_str(&format!(
            "- `{}`: severity={}, direction={}, organelles={}\n",
            signal.id,
            severity_str(signal.severity),
            direction_str(signal.direction),
            signal
                .organelles
                .iter()
                .map(|o| organelle_label(*o))
                .collect::<Vec<_>>()
                .join(", ")
        ));
        for evidence in signal.evidence.iter().take(2) {
            out.push_str(&format!(
                "  - evidence: source={}, organelle={}, axis={}, metric={}, value={:.3}\n",
                source_str(evidence.source),
                organelle_label(evidence.organelle),
                evidence.axis,
                metric_str(evidence.metric),
                evidence.value
            ));
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
    out.push_str(&format!("### {}\n", item.tool));
    match state.organelle_states.iter().find(|s| s.tool == item.tool) {
        Some(organelle) => {
            out.push_str(&format!(
                "- Organelle: `{}`\n",
                organelle_label(organelle.organelle)
            ));
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

    let tool_issues = state
        .issues
        .iter()
        .filter(|i| i.tool.as_deref() == Some(item.tool))
        .collect::<Vec<_>>();
    out.push_str(&format!("- Issues: `{}`\n\n", tool_issues.len()));
}

fn append_axes(out: &mut String, organelle: &OrganelleState) {
    if organelle.axes.is_empty() {
        out.push_str("- Axis medians: not available.\n");
        return;
    }
    let top_axes = organelle
        .axes
        .iter()
        .take(3)
        .map(|axis| format!("{}={:.3}", axis.name, axis.median))
        .collect::<Vec<_>>()
        .join(", ");
    out.push_str(&format!("- Axis medians (first 3): {top_axes}\n"));
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
    if let Some((name, frac)) = dominant {
        out.push_str(&format!("- Dominant regime: {} ({:.3})\n", name, frac));
    } else {
        out.push_str("- Dominant regime: not available.\n");
    }
}

fn append_qc(out: &mut String, organelle: &OrganelleState) {
    if organelle.qc.is_empty() {
        out.push_str("- QC fractions: not available.\n");
        return;
    }
    let qc = organelle
        .qc
        .iter()
        .take(3)
        .map(|(k, v)| format!("{k}={v:.3}"))
        .collect::<Vec<_>>()
        .join(", ");
    out.push_str(&format!("- QC fractions (first 3): {qc}\n"));
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
