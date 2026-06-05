use kira_organelle::report::fii::cai::compute_cai;

fn linspace(n: usize, start: f64, end: f64) -> Vec<f64> {
    if n <= 1 {
        return vec![start];
    }
    (0..n)
        .map(|i| start + (end - start) * (i as f64 / (n as f64 - 1.0)))
        .collect()
}

#[test]
fn synthetic_distributions_show_expected_ordering() {
    let uniform = linspace(1000, 0.2, 0.8);
    let mut skewed = vec![0.05; 900];
    skewed.extend(vec![0.95; 100]);
    let mut bimodal = vec![0.15; 500];
    bimodal.extend(vec![0.85; 500]);

    let a = compute_cai("u", 0, &uniform, 0.66, 0.33, 0.33, 0.34);
    let b = compute_cai("s", 1, &skewed, 0.66, 0.33, 0.33, 0.34);
    let c = compute_cai("b", 2, &bimodal, 0.66, 0.33, 0.33, 0.34);

    assert!(b.cai > a.cai);
    assert!(b.cai > c.cai || c.cai > a.cai);
}

#[test]
fn sensitive_to_small_hard_tail() {
    let diffuse = vec![0.55; 1000];
    let mut small_tail = vec![0.55; 980];
    small_tail.extend(vec![0.95; 20]);

    let a = compute_cai("d", 0, &diffuse, 0.66, 0.33, 0.33, 0.34);
    let b = compute_cai("t", 1, &small_tail, 0.66, 0.33, 0.33, 0.34);
    assert!(b.cai > a.cai);
}

#[test]
fn stable_under_count_scaling() {
    let base = vec![0.2, 0.2, 0.3, 0.4, 0.9];
    let mut scaled = Vec::new();
    for _ in 0..200 {
        scaled.extend(base.iter().copied());
    }
    let a = compute_cai("a", 0, &base, 0.66, 0.33, 0.33, 0.34);
    let b = compute_cai("b", 1, &scaled, 0.66, 0.33, 0.33, 0.34);
    assert!((a.cai - b.cai).abs() < 1e-9);
}

#[test]
fn deterministic_across_runs() {
    let mut values = vec![0.1; 700];
    values.extend(vec![0.95; 40]);
    values.extend(vec![0.5; 260]);
    let a = compute_cai("x", 0, &values, 0.66, 0.33, 0.33, 0.34);
    let b = compute_cai("x", 0, &values, 0.66, 0.33, 0.33, 0.34);
    assert!((a.cai - b.cai).abs() < 1e-12);
    assert!((a.skewness_component - b.skewness_component).abs() < 1e-12);
    assert!((a.tail_heaviness_component - b.tail_heaviness_component).abs() < 1e-12);
    assert!((a.tail_mass - b.tail_mass).abs() < 1e-12);
}
