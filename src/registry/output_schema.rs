pub const BASE_COLUMNS: &[&str] = &["cell_id", "cluster"];

pub const COMPOSITE_COLUMNS: &[&str] = &[
    "StressVector",
    "CompensationDeficit",
    "CPI",
    "AFS",
    "IMSC",
    "RegimeClass",
];

pub const AGG_COLUMNS: &[&str] = &["RareState", "MahalanobisDistance"];

pub const FRAGILITY_COLUMNS: &[&str] = &["DominantAxis", "DominantAxisSensitivity"];
pub const LANDSCAPE_COLUMNS: &[&str] = &[
    "Potential",
    "StabilityGradient",
    "TPI_landscape",
    "BasinId",
    "TransitionCandidate",
];
pub const THERAPEUTIC_COLUMNS: &[&str] = &["DominantVulnerableAxis"];
pub const DYNAMIC_COLUMNS: &[&str] = &["LSI", "MaxTrajectory", "BEE"];

#[derive(Debug, Clone, Copy, Default)]
pub struct OutputSchemaAvailability {
    pub include_aggregate: bool,
    pub include_fragility: bool,
    pub include_landscape: bool,
    pub include_therapeutic: bool,
    pub include_dynamic: bool,
}

pub fn metrics_tsv_columns(availability: OutputSchemaAvailability) -> Vec<&'static str> {
    let mut cols = Vec::with_capacity(
        BASE_COLUMNS.len()
            + COMPOSITE_COLUMNS.len()
            + AGG_COLUMNS.len()
            + FRAGILITY_COLUMNS.len()
            + LANDSCAPE_COLUMNS.len()
            + THERAPEUTIC_COLUMNS.len()
            + DYNAMIC_COLUMNS.len(),
    );
    cols.extend_from_slice(BASE_COLUMNS);
    cols.extend_from_slice(COMPOSITE_COLUMNS);
    if availability.include_aggregate {
        cols.extend_from_slice(AGG_COLUMNS);
    }
    if availability.include_fragility {
        cols.extend_from_slice(FRAGILITY_COLUMNS);
    }
    if availability.include_landscape {
        cols.extend_from_slice(LANDSCAPE_COLUMNS);
    }
    if availability.include_therapeutic {
        cols.extend_from_slice(THERAPEUTIC_COLUMNS);
    }
    if availability.include_dynamic {
        cols.extend_from_slice(DYNAMIC_COLUMNS);
    }
    cols
}
