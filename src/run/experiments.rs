use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub enum ExperimentSource {
    Directory(PathBuf),
    RawCountsFile(PathBuf),
    TenxGroup {
        matrix: PathBuf,
        barcodes: PathBuf,
        features_or_genes: PathBuf,
    },
}

#[derive(Debug, Clone)]
pub struct ExperimentSpec {
    pub name: String,
    pub source: ExperimentSource,
}

impl ExperimentSpec {
    pub fn uses_original_directory(&self) -> bool {
        matches!(self.source, ExperimentSource::Directory(_))
    }
}

pub fn discover_experiments(input: &Path) -> Result<Vec<ExperimentSpec>, String> {
    if input.is_file() {
        let name = experiment_name_from_file_path(input);
        return Ok(vec![ExperimentSpec {
            name,
            source: ExperimentSource::RawCountsFile(input.to_path_buf()),
        }]);
    }

    if !input.is_dir() {
        return Err(format!(
            "input path does not exist or is not directory/file: {}",
            input.display()
        ));
    }

    let mut raw = Vec::new();
    let mut tenx: BTreeMap<String, TenxGroupBuilder> = BTreeMap::new();

    let rd = std::fs::read_dir(input)
        .map_err(|e| format!("failed reading input directory {}: {e}", input.display()))?;
    for entry in rd {
        let entry = entry
            .map_err(|e| format!("failed reading input directory {}: {e}", input.display()))?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(file_name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };

        if is_raw_counts_name(file_name) {
            raw.push(ExperimentSpec {
                name: experiment_name_from_filename(file_name),
                source: ExperimentSource::RawCountsFile(path),
            });
            continue;
        }

        if let Some((prefix, kind)) = parse_tenx_member(file_name) {
            let group = tenx.entry(prefix).or_default();
            match kind {
                TenxMember::Matrix => group.matrix = Some(path),
                TenxMember::Barcodes => group.barcodes = Some(path),
                TenxMember::Features => group.features = Some(path),
                TenxMember::Genes => group.genes = Some(path),
            }
        }
    }

    let mut specs = Vec::new();
    specs.extend(raw);

    for (prefix, group) in tenx {
        let Some(matrix) = group.matrix else { continue };
        let Some(barcodes) = group.barcodes else {
            continue;
        };
        let Some(features_or_genes) = group.features.or(group.genes) else {
            continue;
        };
        specs.push(ExperimentSpec {
            name: sanitize_name(&prefix),
            source: ExperimentSource::TenxGroup {
                matrix,
                barcodes,
                features_or_genes,
            },
        });
    }

    if specs.is_empty() {
        return Ok(vec![ExperimentSpec {
            name: directory_name(input),
            source: ExperimentSource::Directory(input.to_path_buf()),
        }]);
    }

    specs.sort_by(|a, b| a.name.cmp(&b.name));
    ensure_unique_names(&mut specs);
    Ok(specs)
}

pub fn materialize_experiment_input(
    spec: &ExperimentSpec,
    staging_root: &Path,
) -> Result<PathBuf, String> {
    match &spec.source {
        ExperimentSource::Directory(path) => Ok(path.clone()),
        ExperimentSource::RawCountsFile(path) => {
            let exp_dir = recreate_experiment_dir(staging_root, &spec.name)?;
            let target_name = if path
                .extension()
                .and_then(|v| v.to_str())
                .map(|v| v.eq_ignore_ascii_case("gz"))
                .unwrap_or(false)
            {
                "raw_counts.tsv.gz"
            } else {
                "raw_counts.tsv"
            };
            link_or_copy(path, &exp_dir.join(target_name))?;
            Ok(exp_dir)
        }
        ExperimentSource::TenxGroup {
            matrix,
            barcodes,
            features_or_genes,
        } => {
            let exp_dir = recreate_experiment_dir(staging_root, &spec.name)?;

            let matrix_name = canonical_name("matrix.mtx", matrix);
            let barcodes_name = canonical_name("barcodes.tsv", barcodes);
            let features_name = {
                let base = if features_or_genes
                    .file_name()
                    .and_then(|n| n.to_str())
                    .map(|n| n.to_ascii_lowercase().contains("genes.tsv"))
                    .unwrap_or(false)
                {
                    "genes.tsv"
                } else {
                    "features.tsv"
                };
                canonical_name(base, features_or_genes)
            };

            link_or_copy(matrix, &exp_dir.join(matrix_name))?;
            link_or_copy(barcodes, &exp_dir.join(barcodes_name))?;
            link_or_copy(features_or_genes, &exp_dir.join(features_name))?;
            Ok(exp_dir)
        }
    }
}

fn recreate_experiment_dir(staging_root: &Path, name: &str) -> Result<PathBuf, String> {
    std::fs::create_dir_all(staging_root)
        .map_err(|e| format!("failed creating {}: {e}", staging_root.display()))?;
    let exp_dir = staging_root.join(name);
    if exp_dir.exists() {
        std::fs::remove_dir_all(&exp_dir)
            .map_err(|e| format!("failed cleaning {}: {e}", exp_dir.display()))?;
    }
    std::fs::create_dir_all(&exp_dir)
        .map_err(|e| format!("failed creating {}: {e}", exp_dir.display()))?;
    Ok(exp_dir)
}

fn link_or_copy(src: &Path, dst: &Path) -> Result<(), String> {
    if dst.exists() {
        std::fs::remove_file(dst).map_err(|e| format!("failed removing {}: {e}", dst.display()))?;
    }
    match std::fs::hard_link(src, dst) {
        Ok(_) => Ok(()),
        Err(_) => std::fs::copy(src, dst)
            .map_err(|e| format!("failed copying {} -> {}: {e}", src.display(), dst.display()))
            .map(|_| ()),
    }
}

fn canonical_name(base: &str, source: &Path) -> String {
    if source
        .extension()
        .and_then(|v| v.to_str())
        .map(|v| v.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        format!("{base}.gz")
    } else {
        base.to_string()
    }
}

#[derive(Debug, Clone, Copy)]
enum TenxMember {
    Matrix,
    Barcodes,
    Features,
    Genes,
}

#[derive(Debug, Default)]
struct TenxGroupBuilder {
    matrix: Option<PathBuf>,
    barcodes: Option<PathBuf>,
    features: Option<PathBuf>,
    genes: Option<PathBuf>,
}

fn parse_tenx_member(file_name: &str) -> Option<(String, TenxMember)> {
    let lower = file_name.to_ascii_lowercase();
    let patterns = [
        (".matrix.mtx.gz", TenxMember::Matrix),
        (".matrix.mtx", TenxMember::Matrix),
        ("_matrix.mtx.gz", TenxMember::Matrix),
        ("_matrix.mtx", TenxMember::Matrix),
        (".barcodes.tsv.gz", TenxMember::Barcodes),
        (".barcodes.tsv", TenxMember::Barcodes),
        ("_barcodes.tsv.gz", TenxMember::Barcodes),
        ("_barcodes.tsv", TenxMember::Barcodes),
        (".features.tsv.gz", TenxMember::Features),
        (".features.tsv", TenxMember::Features),
        ("_features.tsv.gz", TenxMember::Features),
        ("_features.tsv", TenxMember::Features),
        (".genes.tsv.gz", TenxMember::Genes),
        (".genes.tsv", TenxMember::Genes),
        ("_genes.tsv.gz", TenxMember::Genes),
        ("_genes.tsv", TenxMember::Genes),
    ];

    for (suffix, kind) in patterns {
        if !lower.ends_with(suffix) {
            continue;
        }
        let prefix = &file_name[..file_name.len() - suffix.len()];
        if prefix.is_empty() {
            return None;
        }
        return Some((sanitize_name(prefix), kind));
    }
    None
}

fn is_raw_counts_name(file_name: &str) -> bool {
    let lower = file_name.to_ascii_lowercase();
    lower == "raw_counts.tsv"
        || lower == "raw_counts.tsv.gz"
        || lower.ends_with("_raw_counts.tsv")
        || lower.ends_with("_raw_counts.tsv.gz")
        || lower.ends_with(".raw_counts.tsv")
        || lower.ends_with(".raw_counts.tsv.gz")
}

fn experiment_name_from_file_path(path: &Path) -> String {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("experiment");
    experiment_name_from_filename(name)
}

fn experiment_name_from_filename(file_name: &str) -> String {
    let no_gz = file_name.strip_suffix(".gz").unwrap_or(file_name);
    let no_tsv = no_gz.strip_suffix(".tsv").unwrap_or(no_gz);
    sanitize_name(no_tsv)
}

fn directory_name(path: &Path) -> String {
    let raw = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("experiment");
    sanitize_name(raw)
}

pub fn sanitize_name(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for ch in raw.chars() {
        if ch == '/' || ch == '\\' || ch == ':' || ch == '\0' {
            out.push('_');
        } else {
            out.push(ch);
        }
    }
    let trimmed = out.trim();
    if trimmed.is_empty() {
        "experiment".to_string()
    } else {
        trimmed.to_string()
    }
}

fn ensure_unique_names(specs: &mut [ExperimentSpec]) {
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    for spec in specs {
        let n = seen.entry(spec.name.clone()).or_insert(0);
        if *n == 0 {
            *n = 1;
            continue;
        }
        let suffix = *n;
        *n += 1;
        spec.name = format!("{}_{suffix}", spec.name);
    }
}
