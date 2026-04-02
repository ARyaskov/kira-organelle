pub const ENGINE_NAME: &str = "kira-organelle";
pub const MODEL_VERSION: &str = "2.0.0";
pub const NORMALIZATION_VERSION: &str = "median_mad_v1";
pub const COMPOSITES_VERSION: &str = "systems_v1";
pub const FRAGILITY_VERSION: &str = "fdiff_v1";

pub fn tool_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
