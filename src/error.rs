#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorKind {
    Io,
    Parse,
    ContractViolation,
    ToolFailure,
    Internal,
}

pub fn stable_error_code(kind: ErrorKind) -> &'static str {
    match kind {
        ErrorKind::Io => "E_IO",
        ErrorKind::Parse => "E_PARSE",
        ErrorKind::ContractViolation => "E_CONTRACT",
        ErrorKind::ToolFailure => "E_TOOL",
        ErrorKind::Internal => "E_INTERNAL",
    }
}

pub fn classify_issue_code(code: &str) -> ErrorKind {
    if code.contains("PARSE") {
        return ErrorKind::Parse;
    }
    if code.contains("MISSING") || code.contains("CONTRACT") || code.contains("INVALID") {
        return ErrorKind::ContractViolation;
    }
    if code.contains("TOOL") {
        return ErrorKind::ToolFailure;
    }
    if code.contains("IO") || code.contains("READ") || code.contains("WRITE") {
        return ErrorKind::Io;
    }
    ErrorKind::Internal
}

pub fn exit_code_for_kind(kind: ErrorKind) -> i32 {
    match kind {
        ErrorKind::Io => 2,
        ErrorKind::Parse => 3,
        ErrorKind::ContractViolation => 4,
        ErrorKind::ToolFailure => 5,
        ErrorKind::Internal => 1,
    }
}

pub fn classify_error_message(message: &str) -> ErrorKind {
    let upper = message.to_ascii_uppercase();
    if upper.contains("PARSE") {
        return ErrorKind::Parse;
    }
    if upper.contains("MISSING") || upper.contains("CONTRACT") || upper.contains("INVALID") {
        return ErrorKind::ContractViolation;
    }
    if upper.contains("TOOL") || upper.contains("SPAWN") || upper.contains("EXIT ") {
        return ErrorKind::ToolFailure;
    }
    if upper.contains("IO") || upper.contains("READ") || upper.contains("WRITE") {
        return ErrorKind::Io;
    }
    ErrorKind::Internal
}
