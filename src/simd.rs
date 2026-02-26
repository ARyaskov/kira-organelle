#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimdMode {
    Scalar,
    Avx2,
    Neon,
}

pub fn detect_simd_mode() -> SimdMode {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            return SimdMode::Avx2;
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            return SimdMode::Neon;
        }
    }

    SimdMode::Scalar
}

pub fn log_simd_mode() {
    let mode = detect_simd_mode();
    tracing::info!(mode = ?mode, "runtime SIMD mode selected");
}
