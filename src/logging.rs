use time::format_description::well_known::Rfc3339;
use tracing::Level;
use tracing_subscriber::fmt::time::UtcTime;

use crate::cli::LogLevel;

pub fn init(level: LogLevel, no_color: bool) {
    let parsed = parse_level(level);
    let timer = UtcTime::new(Rfc3339);

    tracing_subscriber::fmt()
        .with_timer(timer)
        .with_target(false)
        .with_ansi(!no_color)
        .with_max_level(parsed)
        .init();
}

fn parse_level(level: LogLevel) -> Level {
    match level {
        LogLevel::Error => Level::ERROR,
        LogLevel::Warn => Level::WARN,
        LogLevel::Info => Level::INFO,
        LogLevel::Debug => Level::DEBUG,
        LogLevel::Trace => Level::TRACE,
    }
}
