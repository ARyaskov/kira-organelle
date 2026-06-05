use std::env;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub enum ToolInvocationMode {
    Library,
    Binary(PathBuf),
    Unavailable,
}

pub fn resolve_tool_invocation(tool: &str) -> ToolInvocationMode {
    if let Some(path) = resolve_binary(tool) {
        return ToolInvocationMode::Binary(path);
    }
    ToolInvocationMode::Unavailable
}

fn resolve_binary(tool: &str) -> Option<PathBuf> {
    let env_key = format!(
        "KIRA_ORGANELLE_BIN_{}",
        tool.replace('-', "_").to_ascii_uppercase()
    );
    if let Some(v) = env::var_os(&env_key) {
        let p = PathBuf::from(v);
        if let Some(resolved) = match_executable(&p) {
            return Some(resolved);
        }
    }

    if let Ok(exe) = env::current_exe()
        && let Some(bin_dir) = exe.parent()
        && let Some(resolved) = match_executable(&bin_dir.join(tool))
    {
        return Some(resolved);
    }

    if let Some(path) = find_in_path(tool) {
        return Some(path);
    }

    find_in_neighbor_repo(tool)
}

fn find_in_path(tool: &str) -> Option<PathBuf> {
    let path_env = env::var_os("PATH")?;
    for dir in env::split_paths(&path_env) {
        if let Some(resolved) = match_executable(&dir.join(tool)) {
            return Some(resolved);
        }
    }
    None
}

/// Resolve `path` to an existing executable file, retrying with the OS exe suffix
/// (`.exe` on Windows, empty elsewhere). Returns the canonical match or `None`.
fn match_executable(path: &Path) -> Option<PathBuf> {
    if is_executable_file(path) {
        return Some(path.to_path_buf());
    }
    let suffix = env::consts::EXE_SUFFIX;
    if !suffix.is_empty() {
        let already_has_suffix = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case(suffix.trim_start_matches('.')))
            .unwrap_or(false);
        if !already_has_suffix {
            let candidate: PathBuf = {
                let mut s = path.as_os_str().to_owned();
                s.push(suffix);
                PathBuf::from(s)
            };
            if is_executable_file(&candidate) {
                return Some(candidate);
            }
        }
    }
    None
}

fn is_executable_file(path: &Path) -> bool {
    path.is_file()
}

fn find_in_neighbor_repo(tool: &str) -> Option<PathBuf> {
    if !tool.starts_with("kira-") {
        return None;
    }
    let roots = collect_neighbor_roots();

    for root in roots {
        let repo_dir = root.join(tool);
        if let Some(resolved) =
            match_executable(&repo_dir.join("target").join("release").join(tool))
        {
            return Some(resolved);
        }
        if let Some(resolved) = match_executable(&repo_dir.join("target").join("debug").join(tool))
        {
            return Some(resolved);
        }
    }

    None
}

fn collect_neighbor_roots() -> Vec<PathBuf> {
    let mut roots: Vec<PathBuf> = Vec::new();

    if let Some(root) = env::var_os("KIRA_PIPELINE_ROOT").map(PathBuf::from) {
        roots.push(root);
    }

    if let Ok(cwd) = env::current_dir() {
        roots.extend(cwd.ancestors().map(Path::to_path_buf));
    }

    if let Ok(exe) = env::current_exe()
        && let Some(bin_dir) = exe.parent()
    {
        roots.extend(bin_dir.ancestors().map(Path::to_path_buf));
    }

    dedup_paths(roots)
}

fn dedup_paths(paths: Vec<PathBuf>) -> Vec<PathBuf> {
    let mut unique: Vec<PathBuf> = Vec::new();
    for path in paths {
        if !unique.iter().any(|p| p == &path) {
            unique.push(path);
        }
    }
    unique
}
