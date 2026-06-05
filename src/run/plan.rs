use std::fmt::Write as _;
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
    let mut out = String::with_capacity(512);
    out.push_str("Execution Plan\n");
    let _ = writeln!(&mut out, "input: {}", plan.input.display());
    let _ = writeln!(&mut out, "out_root: {}", plan.out_root.display());
    let _ = writeln!(&mut out, "organelle_out: {}", plan.organelle_out.display());

    for (idx, step) in plan.steps.iter().enumerate() {
        let _ = writeln!(&mut out, "{}. {}", idx + 1, step.tool);
        let _ = writeln!(&mut out, "   input: {}", step.input.display());
        let _ = writeln!(&mut out, "   out: {}", step.out_dir.display());
        let _ = writeln!(&mut out, "   mode: {}", mode_label(&step.mode));
        out.push_str("   threads: ");
        match step.threads {
            Some(v) => {
                let _ = writeln!(&mut out, "{v}");
            }
            None => out.push_str("-\n"),
        }
        out.push_str("   cache: ");
        match step.cache_path.as_ref() {
            Some(v) => {
                let _ = writeln!(&mut out, "{}", v.display());
            }
            None => out.push_str("disabled\n"),
        }
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
