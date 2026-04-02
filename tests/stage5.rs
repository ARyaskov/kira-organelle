use std::fs;

use kira_organelle::cli::RunArgs;
use kira_organelle::run::exec::{execute_step, prepare_layout};
use kira_organelle::run::plan::{ToolStep, build_execution_plan, render_dry_run_plan};
use kira_organelle::run::tool::ToolInvocationMode;
use tempfile::tempdir;

#[test]
fn dry_run_plan_deterministic() {
    let dir = tempdir().expect("tempdir");
    let args = RunArgs {
        input: dir.path().join("input"),
        integration_manifest: None,
        out: Some(dir.path().join("out")),
        threads: Some(4),
        no_cache: false,
        strict: false,
        dry_run: true,
        fii_weights: None,
        mitoqc_write_profile_json: false,
        write_cells_json: false,
    };

    let plan_a = build_execution_plan(&args);
    let plan_b = build_execution_plan(&args);

    let text_a = render_dry_run_plan(&plan_a);
    let text_b = render_dry_run_plan(&plan_b);

    assert_eq!(text_a, text_b);
}

#[test]
fn binary_fallback_mock() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    let out = dir.path().join("kira-mitoqc");
    fs::create_dir_all(&input).expect("mkdir input");
    fs::create_dir_all(&out).expect("mkdir out");

    let script = dir.path().join("fake-kira-mitoqc.sh");
    fs::write(
        &script,
        "#!/bin/sh\nset -e\nout=''\nwhile [ $# -gt 0 ]; do\n  if [ \"$1\" = \"--out\" ]; then out=\"$2\"; shift 2; continue; fi\n  shift\ndone\ntouch \"$out/mock-ok.txt\"\n",
    )
    .expect("write script");

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = fs::metadata(&script).expect("meta").permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&script, perms).expect("chmod");
    }

    let step = ToolStep {
        tool: "kira-mitoqc".to_string(),
        input,
        out_dir: out.clone(),
        threads: Some(2),
        cache_path: None,
        mode: ToolInvocationMode::Binary(script),
    };

    let result = execute_step(&step);
    assert!(result.success, "binary step failed: {}", result.message);
    assert!(out.join("mock-ok.txt").is_file());
}

#[test]
fn binary_step_passes_assets_to_mitoqc() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    let out = dir.path().join("kira-mitoqc");
    let bin_dir = dir.path().join("bin");
    let assets_dir = bin_dir.join("assets");
    let args_log = dir.path().join("args-mito.txt");
    fs::create_dir_all(&input).expect("mkdir input");
    fs::create_dir_all(&out).expect("mkdir out");
    fs::create_dir_all(&assets_dir).expect("mkdir assets");

    let script = bin_dir.join("fake-kira-mitoqc.sh");
    fs::write(
        &script,
        format!(
            "#!/bin/sh\nset -e\nprintf '%s\\n' \"$@\" > \"{}\"\n",
            args_log.display()
        ),
    )
    .expect("write script");

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = fs::metadata(&script).expect("meta").permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&script, perms).expect("chmod");
    }

    let step = ToolStep {
        tool: "kira-mitoqc".to_string(),
        input,
        out_dir: out,
        threads: Some(2),
        cache_path: None,
        mode: ToolInvocationMode::Binary(script),
    };

    let result = execute_step(&step);
    assert!(result.success, "binary step failed: {}", result.message);

    let args = fs::read_to_string(args_log).expect("read args");
    assert!(
        args.contains("--assets"),
        "missing --assets in args: {args}"
    );
    assert!(
        args.contains(&assets_dir.display().to_string()),
        "missing assets path in args: {args}"
    );
    assert!(args.contains("--redox"), "missing --redox in args: {args}");
}

#[test]
fn binary_step_passes_cache_to_nuclearqc() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    let out = dir.path().join("kira-nuclearqc");
    let cache = dir.path().join("kira-organelle.bin");
    let args_log = dir.path().join("args.txt");
    fs::create_dir_all(&input).expect("mkdir input");
    fs::create_dir_all(&out).expect("mkdir out");

    let script = dir.path().join("fake-kira-nuclearqc.sh");
    fs::write(
        &script,
        format!(
            "#!/bin/sh\nset -e\nprintf '%s\\n' \"$@\" > \"{}\"\n",
            args_log.display()
        ),
    )
    .expect("write script");

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = fs::metadata(&script).expect("meta").permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&script, perms).expect("chmod");
    }

    let step = ToolStep {
        tool: "kira-nuclearqc".to_string(),
        input,
        out_dir: out,
        threads: Some(2),
        cache_path: Some(cache.clone()),
        mode: ToolInvocationMode::Binary(script),
    };

    let result = execute_step(&step);
    assert!(result.success, "binary step failed: {}", result.message);

    let args = fs::read_to_string(args_log).expect("read args");
    assert!(args.contains("--cache"), "missing --cache in args: {args}");
    assert!(
        args.contains(&cache.display().to_string()),
        "missing cache path in args: {args}"
    );
}

#[test]
fn run_layout_created() {
    let dir = tempdir().expect("tempdir");
    let args = RunArgs {
        input: dir.path().join("input"),
        integration_manifest: None,
        out: Some(dir.path().join("out")),
        threads: None,
        no_cache: true,
        strict: false,
        dry_run: false,
        fii_weights: None,
        mitoqc_write_profile_json: false,
        write_cells_json: false,
    };

    let plan = build_execution_plan(&args);
    prepare_layout(&plan).expect("prepare layout");

    assert!(plan.out_root.is_dir());
    assert!(plan.organelle_out.is_dir());
    for step in &plan.steps {
        assert!(step.out_dir.is_dir(), "missing dir for {}", step.tool);
    }
}

#[test]
fn default_cache_path_is_input_root() {
    let dir = tempdir().expect("tempdir");
    let input = dir.path().join("input");
    fs::create_dir_all(&input).expect("mkdir input");
    let args = RunArgs {
        input: input.clone(),
        integration_manifest: None,
        out: Some(dir.path().join("out")),
        threads: None,
        no_cache: false,
        strict: false,
        dry_run: false,
        fii_weights: None,
        mitoqc_write_profile_json: false,
        write_cells_json: false,
    };

    let plan = build_execution_plan(&args);
    for step in &plan.steps {
        assert_eq!(
            step.cache_path.as_deref(),
            Some(input.join("kira-organelle.bin").as_path())
        );
    }
}
