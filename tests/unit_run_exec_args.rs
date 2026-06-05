use std::path::PathBuf;

use kira_organelle::run::exec::build_tool_args;
use kira_organelle::run::plan::ToolStep;
use kira_organelle::run::tool::ToolInvocationMode;
use tempfile::tempdir;

fn build_step(tool: &str, input: &std::path::Path, out_dir: &std::path::Path) -> ToolStep {
    ToolStep {
        tool: tool.to_string(),
        input: input.to_path_buf(),
        out_dir: out_dir.to_path_buf(),
        threads: None,
        cache_path: None,
        mode: ToolInvocationMode::Unavailable,
    }
}

fn out_value(args: &[String]) -> &str {
    let idx = args.iter().position(|v| v == "--out").expect("--out flag");
    args[idx + 1].as_str()
}

#[test]
fn build_tool_args_riboqc_keeps_step_out() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    let out_dir = dir.path().join("out").join("kira-riboqc");
    let bin = dir.path().join("bin").join("kira-riboqc");

    let step = build_step("kira-riboqc", &input, &out_dir);
    let args = build_tool_args(&step, &bin).expect("args");
    assert_eq!(out_value(&args), out_dir.display().to_string());
}

#[test]
fn build_tool_args_riboqc_keeps_root_out_when_already_root() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    let out_dir = dir.path().join("out");
    let bin = dir.path().join("bin").join("kira-riboqc");

    let step = build_step("kira-riboqc", &input, &out_dir);
    let args = build_tool_args(&step, &bin).expect("args");
    assert_eq!(out_value(&args), out_dir.display().to_string());
}

#[test]
fn build_tool_args_non_riboqc_keeps_step_out() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    let out_dir = dir.path().join("out").join("kira-proteoqc");
    let bin = dir.path().join("bin").join("kira-proteoqc");

    let step = build_step("kira-proteoqc", &input, &out_dir);
    let args = build_tool_args(&step, &bin).expect("args");
    assert_eq!(out_value(&args), out_dir.display().to_string());
}

#[test]
fn build_tool_args_input_is_passed_through() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    let out_dir = dir.path().join("out").join("kira-proteoqc");
    let bin: PathBuf = dir.path().join("bin").join("kira-proteoqc");

    let step = build_step("kira-proteoqc", &input, &out_dir);
    let args = build_tool_args(&step, &bin).expect("args");
    let idx = args.iter().position(|v| v == "--input").expect("--input");
    assert!(args[idx + 1].contains("input"));
}
