use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

#[derive(Debug, Parser)]
#[command(name = "kira-organelle")]
#[command(about = "Kira organelle QC aggregator/orchestrator")]
pub struct Cli {
    #[arg(long, global = true, default_value_t = LogLevel::Info, value_enum)]
    pub log_level: LogLevel,
    #[arg(long, global = true)]
    pub no_color: bool,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    Aggregate(AggregateArgs),
    Run(RunArgs),
    ReportFii(ReportFiiArgs),
    ComputeStateDynamics(ComputeStateDynamicsArgs),
    ReportPhase(ReportPhaseArgs),
    Decision(DecisionArgs),
    ReportDecisionTimeline(ReportDecisionTimelineArgs),
    ComputeIli(ComputeIliArgs),
    ComputeCai(ComputeCaiArgs),
    ComputePri(ComputePriArgs),
    ComputeCocs(ComputeCocsArgs),
    ComputeDci(ComputeDciArgs),
}

#[derive(Debug, clap::Args)]
pub struct AggregateArgs {
    #[arg(long)]
    pub input: PathBuf,
    #[arg(long = "input-b")]
    pub input_b: Option<PathBuf>,
    #[arg(long)]
    pub out: Option<PathBuf>,
    #[arg(long)]
    pub strict: bool,
    #[arg(long)]
    pub json: bool,
    #[arg(long = "validate-only")]
    pub validate_only: bool,
    #[arg(long = "fii-weights")]
    pub fii_weights: Option<String>,
}

#[derive(Debug, clap::Args)]
pub struct RunArgs {
    #[arg(long)]
    pub input: PathBuf,
    #[arg(long = "integration-manifest")]
    pub integration_manifest: Option<PathBuf>,
    #[arg(long)]
    pub out: Option<PathBuf>,
    #[arg(long)]
    pub threads: Option<usize>,
    #[arg(long = "no-cache")]
    pub no_cache: bool,
    #[arg(long)]
    pub strict: bool,
    #[arg(long = "dry-run")]
    pub dry_run: bool,
    #[arg(long = "fii-weights")]
    pub fii_weights: Option<String>,
}

#[derive(Debug, clap::Args)]
pub struct ReportFiiArgs {
    #[arg(long, value_delimiter = ',')]
    pub inputs: Vec<PathBuf>,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long)]
    pub manifest: Option<PathBuf>,
    #[arg(long)]
    pub title: Option<String>,
}

#[derive(Debug, clap::Args)]
pub struct ComputeStateDynamicsArgs {
    #[arg(long, value_delimiter = ',')]
    pub inputs: Vec<PathBuf>,
    #[arg(long)]
    pub manifest: Option<PathBuf>,
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, clap::Args)]
pub struct ReportPhaseArgs {
    #[arg(long, value_delimiter = ',')]
    pub inputs: Vec<PathBuf>,
    #[arg(long)]
    pub manifest: Option<PathBuf>,
    #[arg(long)]
    pub config: Option<PathBuf>,
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, clap::Args)]
pub struct DecisionArgs {
    #[arg(long = "inputs")]
    pub input: PathBuf,
    #[arg(long)]
    pub config: Option<PathBuf>,
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, clap::Args)]
pub struct ReportDecisionTimelineArgs {
    #[arg(long = "inputs")]
    pub input: PathBuf,
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, clap::Args)]
pub struct ComputeIliArgs {
    #[arg(long = "inputs")]
    pub input: PathBuf,
    #[arg(long)]
    pub config: Option<PathBuf>,
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, clap::Args)]
pub struct ComputeCaiArgs {
    #[arg(long, value_delimiter = ',')]
    pub inputs: Vec<PathBuf>,
    #[arg(long)]
    pub manifest: Option<PathBuf>,
    #[arg(long)]
    pub config: Option<PathBuf>,
    #[arg(long = "cai-weights")]
    pub cai_weights: Option<String>,
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, clap::Args)]
pub struct ComputePriArgs {
    #[arg(long, value_delimiter = ',')]
    pub inputs: Vec<PathBuf>,
    #[arg(long)]
    pub manifest: Option<PathBuf>,
    #[arg(long)]
    pub config: Option<PathBuf>,
    #[arg(long = "pri-weights")]
    pub pri_weights: Option<String>,
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, clap::Args)]
pub struct ComputeCocsArgs {
    #[arg(long = "inputs")]
    pub input: PathBuf,
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, clap::Args)]
pub struct ComputeDciArgs {
    #[arg(long = "inputs")]
    pub input: PathBuf,
    #[arg(long)]
    pub config: Option<PathBuf>,
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}
