pub fn median_in_place(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let n = values.len();
    if n % 2 == 1 {
        let mid = n / 2;
        let (_, pivot, _) = values.select_nth_unstable_by(mid, |a, b| a.total_cmp(b));
        return Some(*pivot);
    }

    let hi = n / 2;
    let (_, pivot_hi, _) = values.select_nth_unstable_by(hi, |a, b| a.total_cmp(b));
    let hi_val = *pivot_hi;
    let mut lo_val = values[0];
    for v in &values[..hi] {
        if v.total_cmp(&lo_val).is_gt() {
            lo_val = *v;
        }
    }
    Some((lo_val + hi_val) / 2.0)
}

pub fn percentile_nearest_rank_in_place(values: &mut [f64], p: f64) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let rank = ((values.len() as f64) * p).ceil() as usize;
    let idx = rank.saturating_sub(1).min(values.len().saturating_sub(1));
    let (_, pivot, _) = values.select_nth_unstable_by(idx, |a, b| a.total_cmp(b));
    Some(*pivot)
}
