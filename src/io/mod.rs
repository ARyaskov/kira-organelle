use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use serde_json::Value;

pub const DETERMINISTIC_TIMESTAMP_RFC3339: &str = "1970-01-01T00:00:00Z";

pub fn now_rfc3339_utc() -> String {
    let now = time::OffsetDateTime::now_utc();
    now.format(&time::format_description::well_known::Rfc3339)
        .unwrap_or_else(|_| DETERMINISTIC_TIMESTAMP_RFC3339.to_string())
}

pub fn deterministic_rfc3339_utc() -> String {
    DETERMINISTIC_TIMESTAMP_RFC3339.to_string()
}

pub fn write_json_atomic(path: &Path, value: &Value) -> Result<(), std::io::Error> {
    let mut bytes = serde_json::to_vec_pretty(value)
        .map_err(|e| std::io::Error::other(format!("json serialize error: {e}")))?;
    bytes.push(b'\n');
    write_bytes_atomic(path, &bytes)
}

pub fn write_bytes_atomic(path: &Path, bytes: &[u8]) -> Result<(), std::io::Error> {
    let parent = path
        .parent()
        .ok_or_else(|| std::io::Error::other("no parent directory"))?;
    fs::create_dir_all(parent)?;

    let tmp_path = parent.join(format!(
        ".{}.tmp",
        path.file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("tmp.out")
    ));

    let mut tmp_file = File::create(&tmp_path)?;
    tmp_file.write_all(bytes)?;
    maybe_sync_file(&tmp_file)?;
    drop(tmp_file);

    fs::rename(&tmp_path, path)?;
    maybe_sync_dir(parent)?;
    Ok(())
}

fn maybe_sync_file(_file: &File) -> Result<(), std::io::Error> {
    #[cfg(feature = "strict-fsync")]
    {
        _file.sync_all()?;
    }
    Ok(())
}

fn maybe_sync_dir(_path: &Path) -> Result<(), std::io::Error> {
    #[cfg(feature = "strict-fsync")]
    {
        let dir = File::open(_path)?;
        dir.sync_all()?;
    }
    Ok(())
}
