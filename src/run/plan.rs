use std::path::{Path, PathBuf};

use crate::cli::RunArgs;

use super::tool::{ToolInvocationMode, resolve_tool_invocation};

pub const PIPELINE_ORDER: [&str; 9] = [
    "kira-mitoqc",
    "kira-nuclearqc",
    "kira-spliceqc",
    "kira-riboqc",
    "kira-proteoqc",
    "kira-autolys",
    "kira-secretion",
    "kira-energetics",
    "kira-microenvironment",
];

#[derive(Debug, Clone)]
pub struct ToolStep {
    pub tool: String,
    pub input: PathBuf,
    pub out_dir: PathBuf,
    pub threads: Option<usize>,
    pub cache_path: Option<PathBuf>,
    pub mode: ToolInvocationMode,
}

#[derive(Debug, Clone)]
pub struct ExecutionPlan {
    pub input: PathBuf,
    pub out_root: PathBuf,
    pub organelle_out: PathBuf,
    pub steps: Vec<ToolStep>,
}

pub fn build_execution_plan(args: &RunArgs) -> ExecutionPlan {
    let out_root = args
        .out
        .clone()
        .unwrap_or_else(|| args.input.join("kira-pipeline"));
    let organelle_out = out_root.join("kira-organelle");

    let cache_path = if args.no_cache {
        None
    } else {
        let input_root = if args.input.is_file() {
            args.input.parent().unwrap_or(&args.input)
        } else {
            &args.input
        };
        Some(input_root.join("kira-organelle.bin"))
    };

    let mut steps = Vec::new();
    for tool in PIPELINE_ORDER {
        steps.push(ToolStep {
            tool: tool.to_string(),
            input: args.input.clone(),
            out_dir: out_root.join(tool),
            threads: args.threads,
            cache_path: cache_path.clone(),
            mode: resolve_tool_invocation(tool),
        });
    }

    ExecutionPlan {
        input: args.input.clone(),
        out_root,
        organelle_out,
        steps,
    }
}

pub fn render_dry_run_plan(plan: &ExecutionPlan) -> String {
    let mut out = String::new();
    out.push_str("Execution Plan\n");
    out.push_str(&format!("input: {}\n", plan.input.display()));
    out.push_str(&format!("out_root: {}\n", plan.out_root.display()));
    out.push_str(&format!(
        "organelle_out: {}\n",
        plan.organelle_out.display()
    ));

    for (idx, step) in plan.steps.iter().enumerate() {
        out.push_str(&format!("{}. {}\n", idx + 1, step.tool));
        out.push_str(&format!("   input: {}\n", step.input.display()));
        out.push_str(&format!("   out: {}\n", step.out_dir.display()));
        out.push_str(&format!("   mode: {}\n", mode_label(&step.mode)));
        out.push_str(&format!(
            "   threads: {}\n",
            step.threads
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string())
        ));
        out.push_str(&format!(
            "   cache: {}\n",
            step.cache_path
                .as_ref()
                .map(|v| v.display().to_string())
                .unwrap_or_else(|| "disabled".to_string())
        ));
    }

    out
}

fn mode_label(mode: &ToolInvocationMode) -> String {
    match mode {
        ToolInvocationMode::Library => "library".to_string(),
        ToolInvocationMode::Binary(path) => format!("binary:{}", path.display()),
        ToolInvocationMode::Unavailable => "unavailable".to_string(),
    }
}

#[allow(dead_code)]
fn _assert_path(_: &Path) {}
