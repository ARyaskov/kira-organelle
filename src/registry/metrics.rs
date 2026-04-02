use super::composites::CompositeId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MetricKind {
    RawScore,
    ZScore,
    Flag,
    Category,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ToolId {
    NuclearQc,
    SpliceQc,
    MitoQc,
    Energetics,
    RiboQc,
    ProteoQc,
    Autolys,
    Microenvironment,
}

impl ToolId {
    pub fn from_tool_name(name: &str) -> Option<Self> {
        match name {
            "kira-nuclearqc" => Some(Self::NuclearQc),
            "kira-spliceqc" => Some(Self::SpliceQc),
            "kira-mitoqc" => Some(Self::MitoQc),
            "kira-energetics" => Some(Self::Energetics),
            "kira-riboqc" => Some(Self::RiboQc),
            "kira-proteoqc" => Some(Self::ProteoQc),
            "kira-autolys" => Some(Self::Autolys),
            "kira-microenvironment" => Some(Self::Microenvironment),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(usize)]
pub enum MetricId {
    Rss = 0,
    Ddr = 1,
    Cds = 2,
    Sas = 3,
    Sos = 4,
    Rlr = 5,
    Sii = 6,
    Mri = 7,
    Osl = 8,
    Ess = 9,
    Mcb = 10,
    Ogi = 11,
    Apb = 12,
    Ebs = 13,
    Mfi = 14,
    Lep = 15,
    Imsi = 16,
    Tpi = 17,
    Rbl = 18,
    MtorP = 19,
    IsrA = 20,
    Tpib = 21,
    Tsm = 22,
    Cci = 23,
    Pci = 24,
    UprA = 25,
    Pls = 26,
    Sci = 27,
    Pcp = 28,
    Ais = 29,
    Afs = 30,
    Lci = 31,
    Afp = 32,
    CdsAutolys = 33,
    Asm = 34,
    MitophagyScore = 35,
    Hsi = 36,
    Ias = 37,
    Iss = 38,
    Mio = 39,
    SiiStromal = 40,
    Msm = 41,
}

impl MetricId {
    pub const COUNT: usize = 42;

    pub fn as_index(self) -> usize {
        self as usize
    }
}

#[derive(Debug, Clone, Copy)]
pub struct MetricSpec {
    pub id: MetricId,
    pub canonical_name: &'static str,
    pub accepted_aliases: &'static [&'static str],
    pub source_tool: ToolId,
    pub kind: MetricKind,
    pub required_for: &'static [CompositeId],
}

const REQ_STRESS_VECTOR: &[CompositeId] = &[CompositeId::StressVector, CompositeId::Cpi];
const REQ_COMP_DEF: &[CompositeId] = &[CompositeId::CompensationDeficit, CompositeId::Cpi];
const REQ_AFS: &[CompositeId] = &[CompositeId::Afs];
const REQ_IMSC: &[CompositeId] = &[CompositeId::Imsc];
const REQ_NONE: &[CompositeId] = &[];

pub const METRIC_SPECS: &[MetricSpec] = &[
    MetricSpec {
        id: MetricId::Rss,
        canonical_name: "RSS",
        accepted_aliases: &["rss", "replication_stress_score"],
        source_tool: ToolId::NuclearQc,
        kind: MetricKind::RawScore,
        required_for: REQ_COMP_DEF,
    },
    MetricSpec {
        id: MetricId::Ddr,
        canonical_name: "DDR",
        accepted_aliases: &["ddr", "dna_damage_response"],
        source_tool: ToolId::NuclearQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Cds,
        canonical_name: "CDS",
        accepted_aliases: &["cds", "chromatin_destabilization_score"],
        source_tool: ToolId::NuclearQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Sas,
        canonical_name: "SAS",
        accepted_aliases: &["sas", "senescence_associated_signal"],
        source_tool: ToolId::NuclearQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Sos,
        canonical_name: "SOS",
        accepted_aliases: &["sos", "splicing_overload_score"],
        source_tool: ToolId::SpliceQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Rlr,
        canonical_name: "RLR",
        accepted_aliases: &["rlr", "retention_load_ratio"],
        source_tool: ToolId::SpliceQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Sii,
        canonical_name: "SII",
        accepted_aliases: &[
            "sii",
            "sis",
            "splice_stress_index",
            "splicing_stability_index",
        ],
        source_tool: ToolId::SpliceQc,
        kind: MetricKind::RawScore,
        required_for: REQ_COMP_DEF,
    },
    MetricSpec {
        id: MetricId::Mri,
        canonical_name: "MRI",
        accepted_aliases: &["mri", "mitochondrial_rigidity_index"],
        source_tool: ToolId::MitoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Osl,
        canonical_name: "OSL",
        accepted_aliases: &["osl", "oxidative_stress_load"],
        source_tool: ToolId::MitoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_COMP_DEF,
    },
    MetricSpec {
        id: MetricId::Ess,
        canonical_name: "ESS",
        accepted_aliases: &["ess", "energetic_stall_score"],
        source_tool: ToolId::MitoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Mcb,
        canonical_name: "MCB",
        accepted_aliases: &["mcb", "mitochondrial_compensation_buffer"],
        source_tool: ToolId::MitoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_AFS,
    },
    MetricSpec {
        id: MetricId::Ogi,
        canonical_name: "OGI",
        accepted_aliases: &["ogi", "organelle_gating_instability"],
        source_tool: ToolId::MitoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Apb,
        canonical_name: "APB",
        accepted_aliases: &["apb", "atp_buffer"],
        source_tool: ToolId::Energetics,
        kind: MetricKind::RawScore,
        required_for: REQ_COMP_DEF,
    },
    MetricSpec {
        id: MetricId::Ebs,
        canonical_name: "EBS",
        accepted_aliases: &["ebs", "energetic_barrier_score"],
        source_tool: ToolId::Energetics,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Mfi,
        canonical_name: "MFI",
        accepted_aliases: &["mfi", "metabolic_flexibility_index"],
        source_tool: ToolId::Energetics,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Lep,
        canonical_name: "LEP",
        accepted_aliases: &["lep", "lactate_exhaustion_proxy"],
        source_tool: ToolId::Energetics,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Imsi,
        canonical_name: "IMSI",
        accepted_aliases: &["imsi", "immune_metabolic_suppression_index"],
        source_tool: ToolId::Energetics,
        kind: MetricKind::RawScore,
        required_for: REQ_IMSC,
    },
    MetricSpec {
        id: MetricId::Tpi,
        canonical_name: "TPI",
        accepted_aliases: &["tpi", "translation_pressure_index"],
        source_tool: ToolId::RiboQc,
        kind: MetricKind::RawScore,
        required_for: REQ_COMP_DEF,
    },
    MetricSpec {
        id: MetricId::Rbl,
        canonical_name: "RBL",
        accepted_aliases: &["rbl", "ribosome_burden_load"],
        source_tool: ToolId::RiboQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::MtorP,
        canonical_name: "mTOR_P",
        accepted_aliases: &["mtor_p", "mtor_pressure"],
        source_tool: ToolId::RiboQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::IsrA,
        canonical_name: "ISR_A",
        accepted_aliases: &["isr_a", "integrated_stress_response_activity"],
        source_tool: ToolId::RiboQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Tpib,
        canonical_name: "TPIB",
        accepted_aliases: &["tpib", "translation_pressure_imbalance"],
        source_tool: ToolId::RiboQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Tsm,
        canonical_name: "TSM",
        accepted_aliases: &["tsm", "translation_stress_mean", "translation_stress_index"],
        source_tool: ToolId::RiboQc,
        kind: MetricKind::RawScore,
        required_for: REQ_STRESS_VECTOR,
    },
    MetricSpec {
        id: MetricId::Cci,
        canonical_name: "CCI",
        accepted_aliases: &["cci", "chaperone_capacity_index"],
        source_tool: ToolId::ProteoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_AFS,
    },
    MetricSpec {
        id: MetricId::Pci,
        canonical_name: "PCI",
        accepted_aliases: &["pci", "proteasome_capacity_index"],
        source_tool: ToolId::ProteoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_AFS,
    },
    MetricSpec {
        id: MetricId::UprA,
        canonical_name: "UPR_A",
        accepted_aliases: &["upr_a", "upr_activity"],
        source_tool: ToolId::ProteoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Pls,
        canonical_name: "PLS",
        accepted_aliases: &["pls", "proteostasis_load_score"],
        source_tool: ToolId::ProteoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Sci,
        canonical_name: "SCI",
        accepted_aliases: &["sci", "stress_chaperone_induction"],
        source_tool: ToolId::ProteoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Pcp,
        canonical_name: "PCP",
        accepted_aliases: &["pcp", "proteostasis_collapse_proxy"],
        source_tool: ToolId::ProteoQc,
        kind: MetricKind::RawScore,
        required_for: REQ_STRESS_VECTOR,
    },
    MetricSpec {
        id: MetricId::Ais,
        canonical_name: "AIS",
        accepted_aliases: &["ais", "autophagy_initiation_score"],
        source_tool: ToolId::Autolys,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Afs,
        canonical_name: "AFS",
        accepted_aliases: &["afs", "autophagy_flux_score"],
        source_tool: ToolId::Autolys,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Lci,
        canonical_name: "LCI",
        accepted_aliases: &["lci", "lysosomal_capacity_index"],
        source_tool: ToolId::Autolys,
        kind: MetricKind::RawScore,
        required_for: REQ_AFS,
    },
    MetricSpec {
        id: MetricId::Afp,
        canonical_name: "AFP",
        accepted_aliases: &["afp", "autophagy_flux_proxy"],
        source_tool: ToolId::Autolys,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::CdsAutolys,
        canonical_name: "CDS_autolys",
        accepted_aliases: &["cds_autolys", "cargo_degradation_score"],
        source_tool: ToolId::Autolys,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Asm,
        canonical_name: "ASM",
        accepted_aliases: &["asm", "autophagy_stress_mean", "stress_autophagy_index"],
        source_tool: ToolId::Autolys,
        kind: MetricKind::RawScore,
        required_for: REQ_STRESS_VECTOR,
    },
    MetricSpec {
        id: MetricId::MitophagyScore,
        canonical_name: "MitophagyScore",
        accepted_aliases: &["mitophagyscore", "mitophagy_score"],
        source_tool: ToolId::Autolys,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Hsi,
        canonical_name: "HSI",
        accepted_aliases: &["hsi", "hypoxia_stress_index"],
        source_tool: ToolId::Microenvironment,
        kind: MetricKind::RawScore,
        required_for: REQ_COMP_DEF,
    },
    MetricSpec {
        id: MetricId::Ias,
        canonical_name: "IAS",
        accepted_aliases: &["ias", "immune_activation_score"],
        source_tool: ToolId::Microenvironment,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Iss,
        canonical_name: "ISS",
        accepted_aliases: &["iss", "immune_suppression_score"],
        source_tool: ToolId::Microenvironment,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Mio,
        canonical_name: "MIO",
        accepted_aliases: &["mio", "matrix_interaction_overload"],
        source_tool: ToolId::Microenvironment,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::SiiStromal,
        canonical_name: "SII_stromal",
        accepted_aliases: &["sii_stromal", "stromal_stress_index"],
        source_tool: ToolId::Microenvironment,
        kind: MetricKind::RawScore,
        required_for: REQ_NONE,
    },
    MetricSpec {
        id: MetricId::Msm,
        canonical_name: "MSM",
        accepted_aliases: &["msm", "microenvironment_suppression_metric"],
        source_tool: ToolId::Microenvironment,
        kind: MetricKind::RawScore,
        required_for: &[
            CompositeId::StressVector,
            CompositeId::Cpi,
            CompositeId::Imsc,
            CompositeId::RegimeClass,
        ],
    },
];

pub fn metric_specs() -> &'static [MetricSpec] {
    METRIC_SPECS
}

pub fn find_metric_spec(tool: ToolId, column_name: &str) -> Option<&'static MetricSpec> {
    let exact = METRIC_SPECS
        .iter()
        .find(|spec| spec.source_tool == tool && spec.canonical_name == column_name);
    if exact.is_some() {
        return exact;
    }
    let lower = column_name.to_ascii_lowercase();
    METRIC_SPECS.iter().find(|spec| {
        if spec.source_tool != tool {
            return false;
        }
        spec.accepted_aliases
            .iter()
            .any(|alias| lower.eq_ignore_ascii_case(alias))
    })
}

pub fn expected_metrics_for_tool(tool: ToolId) -> Vec<&'static MetricSpec> {
    METRIC_SPECS
        .iter()
        .filter(|spec| spec.source_tool == tool)
        .collect()
}

pub fn metric_spec_by_id(id: MetricId) -> &'static MetricSpec {
    &METRIC_SPECS[id.as_index()]
}

#[cfg(test)]
mod tests {
    use super::{MetricId, ToolId, find_metric_spec, metric_spec_by_id};

    #[test]
    fn canonical_order_and_count_are_stable() {
        assert_eq!(MetricId::COUNT, 42);
        assert_eq!(metric_spec_by_id(MetricId::Rss).canonical_name, "RSS");
        assert_eq!(metric_spec_by_id(MetricId::Msm).canonical_name, "MSM");
    }

    #[test]
    fn alias_lookup_is_supported() {
        let spec = find_metric_spec(ToolId::RiboQc, "translation_stress_index").expect("alias");
        assert_eq!(spec.canonical_name, "TSM");
        let exact = find_metric_spec(ToolId::RiboQc, "mTOR_P").expect("exact");
        assert_eq!(exact.canonical_name, "mTOR_P");
    }
}
