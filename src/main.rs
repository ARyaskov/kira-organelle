use clap::Parser;
use kira_organelle::cli::{Cli, Commands};

fn main() {
    let cli = Cli::parse();
    kira_organelle::logging::init(cli.log_level, cli.no_color);
    kira_organelle::simd::log_simd_mode();

    match cli.command {
        Commands::Aggregate(args) => {
            let fii_weights = match args.fii_weights.as_deref() {
                Some(raw) => match kira_organelle::fii::FiiWeights::parse(raw) {
                    Ok(v) => Some(v),
                    Err(err) => {
                        tracing::error!("invalid --fii-weights: {err}");
                        std::process::exit(2);
                    }
                },
                None => None,
            };
            let opts = kira_organelle::AggregateOptions {
                input: args.input,
                input_b: args.input_b,
                out: args.out,
                strict: args.strict,
                json: args.json,
                validate_only: args.validate_only,
                fii_weights,
            };
            if let Err(err) = kira_organelle::run_aggregate(&opts) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::Run(args) => {
            if let Err(err) = kira_organelle::run_run_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::ReportFii(args) => {
            if let Err(err) = kira_organelle::run_report_fii_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::ComputeStateDynamics(args) => {
            if let Err(err) = kira_organelle::run_compute_state_dynamics_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::ReportPhase(args) => {
            if let Err(err) = kira_organelle::run_report_phase_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::Decision(args) => {
            if let Err(err) = kira_organelle::run_decision_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::ReportDecisionTimeline(args) => {
            if let Err(err) = kira_organelle::run_report_decision_timeline_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::ComputeIli(args) => {
            if let Err(err) = kira_organelle::run_compute_ili_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::ComputeCai(args) => {
            if let Err(err) = kira_organelle::run_compute_cai_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::ComputePri(args) => {
            if let Err(err) = kira_organelle::run_compute_pri_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::ComputeCocs(args) => {
            if let Err(err) = kira_organelle::run_compute_cocs_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
        Commands::ComputeDci(args) => {
            if let Err(err) = kira_organelle::run_compute_dci_command(&args) {
                let kind = kira_organelle::error::classify_error_message(&err);
                let code = kira_organelle::error::stable_error_code(kind);
                tracing::error!(error_code = %code, "{err}");
                std::process::exit(kira_organelle::error::exit_code_for_kind(kind));
            }
        }
    }
}
