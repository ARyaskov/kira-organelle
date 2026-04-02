use super::metrics::MetricId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CompositeId {
    StressVector,
    CompensationDeficit,
    Cpi,
    Afs,
    Imsc,
    RegimeClass,
    HeterogeneityIndex,
    TailRiskIndex,
    RareStateFraction,
}

#[derive(Debug, Clone, Copy)]
pub struct CompositeSpec {
    pub id: CompositeId,
    pub required_metrics: &'static [MetricId],
    pub optional_metrics: &'static [MetricId],
    pub output_columns: &'static [&'static str],
}

pub const COMPOSITE_SPECS: &[CompositeSpec] = &[
    CompositeSpec {
        id: CompositeId::StressVector,
        required_metrics: &[
            MetricId::Rss,
            MetricId::Sii,
            MetricId::Osl,
            MetricId::Tsm,
            MetricId::Pcp,
            MetricId::Asm,
            MetricId::Msm,
        ],
        optional_metrics: &[],
        output_columns: &["StressVector"],
    },
    CompositeSpec {
        id: CompositeId::CompensationDeficit,
        required_metrics: &[
            MetricId::Tpi,
            MetricId::Cci,
            MetricId::Pci,
            MetricId::Osl,
            MetricId::Lci,
            MetricId::Rss,
            MetricId::Sii,
            MetricId::Hsi,
            MetricId::Apb,
        ],
        optional_metrics: &[],
        output_columns: &["CompensationDeficit"],
    },
    CompositeSpec {
        id: CompositeId::Cpi,
        required_metrics: &[
            MetricId::Rss,
            MetricId::Sii,
            MetricId::Osl,
            MetricId::Tsm,
            MetricId::Pcp,
            MetricId::Asm,
            MetricId::Msm,
            MetricId::Tpi,
            MetricId::Cci,
            MetricId::Pci,
            MetricId::Lci,
            MetricId::Hsi,
            MetricId::Apb,
        ],
        optional_metrics: &[],
        output_columns: &["CPI"],
    },
    CompositeSpec {
        id: CompositeId::Afs,
        required_metrics: &[MetricId::Cci, MetricId::Pci, MetricId::Lci, MetricId::Mcb],
        optional_metrics: &[],
        output_columns: &["AFS"],
    },
    CompositeSpec {
        id: CompositeId::Imsc,
        required_metrics: &[MetricId::Imsi, MetricId::Msm],
        optional_metrics: &[],
        output_columns: &["IMSC"],
    },
    CompositeSpec {
        id: CompositeId::RegimeClass,
        required_metrics: &[],
        optional_metrics: &[],
        output_columns: &["RegimeClass"],
    },
    CompositeSpec {
        id: CompositeId::HeterogeneityIndex,
        required_metrics: &[],
        optional_metrics: &[],
        output_columns: &["HeterogeneityIndex"],
    },
    CompositeSpec {
        id: CompositeId::TailRiskIndex,
        required_metrics: &[],
        optional_metrics: &[],
        output_columns: &["TailRiskIndex"],
    },
    CompositeSpec {
        id: CompositeId::RareStateFraction,
        required_metrics: &[],
        optional_metrics: &[],
        output_columns: &["RareStateFraction"],
    },
];

pub fn composite_specs() -> &'static [CompositeSpec] {
    COMPOSITE_SPECS
}
