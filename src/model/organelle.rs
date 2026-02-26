use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OrganelleId {
    Mitochondria,
    Nucleus,
    Spliceosome,
    Ribosome,
    Proteostasis,
    Autophagy,
    Secretion,
    Energetics,
    Microenvironment,
}

impl OrganelleId {
    pub fn from_tool_name(tool: &str) -> Option<Self> {
        match tool {
            "kira-mitoqc" => Some(Self::Mitochondria),
            "kira-nuclearqc" => Some(Self::Nucleus),
            "kira-spliceqc" => Some(Self::Spliceosome),
            "kira-riboqc" => Some(Self::Ribosome),
            "kira-proteoqc" => Some(Self::Proteostasis),
            "kira-autolys" => Some(Self::Autophagy),
            "kira-secretion" => Some(Self::Secretion),
            "kira-energetics" => Some(Self::Energetics),
            "kira-microenvironment" => Some(Self::Microenvironment),
            _ => None,
        }
    }
}
