use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde_json::{Map, Value, json};

const ABS_DELTA_THRESHOLD: f64 = 0.10;

pub fn write_timeseries_report(
    experiments: &[(String, PathBuf)],
    report_dir: &Path,
) -> Result<(), String> {
    if experiments.len() <= 1 {
        return Ok(());
    }
    std::fs::create_dir_all(report_dir)
        .map_err(|e| format!("failed creating {}: {e}", report_dir.display()))?;

    let data = build_data_json(experiments)?;
    let mut data_bytes = serde_json::to_vec_pretty(&data)
        .map_err(|e| format!("failed serializing timeseries data: {e}"))?;
    data_bytes.push(b'\n');
    crate::io::write_bytes_atomic(&report_dir.join("data.json"), &data_bytes)
        .map_err(|e| format!("failed writing data.json: {e}"))?;

    crate::io::write_bytes_atomic(&report_dir.join("index.html"), INDEX_HTML.as_bytes())
        .map_err(|e| format!("failed writing index.html: {e}"))?;
    crate::io::write_bytes_atomic(&report_dir.join("styles.css"), STYLES_CSS.as_bytes())
        .map_err(|e| format!("failed writing styles.css: {e}"))?;
    crate::io::write_bytes_atomic(&report_dir.join("app.js"), APP_JS.as_bytes())
        .map_err(|e| format!("failed writing app.js: {e}"))?;
    Ok(())
}

fn build_data_json(experiments: &[(String, PathBuf)]) -> Result<Value, String> {
    let mut metric_series: BTreeMap<String, Vec<Option<f64>>> = BTreeMap::new();
    let mut metric_labels: BTreeMap<String, String> = BTreeMap::new();
    let mut significant = Vec::new();

    for (idx, (_name, out_dir)) in experiments.iter().enumerate() {
        let state_path = out_dir.join("kira-organelle").join("state.json");
        let raw = std::fs::read_to_string(&state_path)
            .map_err(|e| format!("failed reading {}: {e}", state_path.display()))?;
        let state: Value = serde_json::from_str(&raw)
            .map_err(|e| format!("failed parsing {}: {e}", state_path.display()))?;
        let metrics = flatten_state_metrics(&state);

        for vals in metric_series.values_mut() {
            while vals.len() < idx {
                vals.push(None);
            }
            vals.push(None);
        }

        for (key, item) in &metrics {
            let series = metric_series
                .entry(key.clone())
                .or_insert_with(|| vec![None; idx + 1]);
            if series.len() < idx + 1 {
                series.resize(idx + 1, None);
            }
            series[idx] = Some(item.value);
            metric_labels
                .entry(key.clone())
                .or_insert_with(|| item.label.clone());
        }
    }

    for vals in metric_series.values_mut() {
        vals.resize(experiments.len(), None);
    }

    for (key, vals) in &metric_series {
        for i in 1..vals.len() {
            let Some(prev) = vals[i - 1] else { continue };
            let Some(curr) = vals[i] else { continue };
            let delta = curr - prev;
            if delta.abs() < ABS_DELTA_THRESHOLD {
                continue;
            }
            significant.push(json!({
                "experiment_index": i,
                "experiment": experiments[i].0,
                "metric_key": key,
                "metric_label": metric_labels.get(key).cloned().unwrap_or_else(|| key.clone()),
                "previous": prev,
                "current": curr,
                "delta": delta,
            }));
        }
    }

    significant.sort_by(|a, b| {
        let ai = a
            .get("experiment_index")
            .and_then(Value::as_u64)
            .unwrap_or(0);
        let bi = b
            .get("experiment_index")
            .and_then(Value::as_u64)
            .unwrap_or(0);
        if ai != bi {
            return ai.cmp(&bi);
        }
        let ad = a.get("delta").and_then(Value::as_f64).unwrap_or(0.0).abs();
        let bd = b.get("delta").and_then(Value::as_f64).unwrap_or(0.0).abs();
        bd.partial_cmp(&ad).unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut metrics_obj = Map::new();
    for (key, vals) in metric_series {
        let label = metric_labels
            .get(&key)
            .cloned()
            .unwrap_or_else(|| key.clone());
        let values = vals
            .into_iter()
            .map(|v| match v {
                Some(n) => Value::from(n),
                None => Value::Null,
            })
            .collect::<Vec<_>>();
        metrics_obj.insert(
            key.clone(),
            json!({
                "label": label,
                "values": values
            }),
        );
    }

    Ok(json!({
        "schema": "kira-organelle-timeseries-v1",
        "experiments": experiments.iter().map(|(name, _)| name.clone()).collect::<Vec<_>>(),
        "significant": significant,
        "metrics": Value::Object(metrics_obj),
    }))
}

#[derive(Debug, Clone)]
struct MetricPoint {
    label: String,
    value: f64,
}

fn flatten_state_metrics(state: &Value) -> BTreeMap<String, MetricPoint> {
    let mut out = BTreeMap::new();
    let Some(orgs) = state.get("organelle_states").and_then(Value::as_array) else {
        return out;
    };
    for org in orgs {
        let organelle = org
            .get("organelle")
            .and_then(Value::as_str)
            .unwrap_or("unknown");
        let tool = org.get("tool").and_then(Value::as_str).unwrap_or(organelle);

        if let Some(axes) = org.get("axes").and_then(Value::as_array) {
            for axis in axes {
                let Some(name) = axis.get("name").and_then(Value::as_str) else {
                    continue;
                };
                let Some(median) = axis.get("median").and_then(Value::as_f64) else {
                    continue;
                };
                let key = format!("{organelle}::axis::{name}");
                out.insert(
                    key,
                    MetricPoint {
                        label: format!("{tool} / axis / {name}"),
                        value: median,
                    },
                );
            }
        }

        if let Some(qc) = org.get("qc").and_then(Value::as_object) {
            for (name, val) in qc {
                let Some(v) = val.as_f64() else { continue };
                let key = format!("{organelle}::qc::{name}");
                out.insert(
                    key,
                    MetricPoint {
                        label: format!("{tool} / qc / {name}"),
                        value: v,
                    },
                );
            }
        }

        if let Some(frac) = org
            .get("regimes")
            .and_then(Value::as_object)
            .and_then(|r| r.get("fractions"))
            .and_then(Value::as_object)
        {
            for (name, val) in frac {
                let Some(v) = val.as_f64() else { continue };
                let key = format!("{organelle}::regime::{name}");
                out.insert(
                    key,
                    MetricPoint {
                        label: format!("{tool} / regime / {name}"),
                        value: v,
                    },
                );
            }
        }
    }
    out
}

const INDEX_HTML: &str = r#"<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Kira Timeseries Report</title>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css">
  <link rel="stylesheet" href="styles.css">
</head>
<body>
  <div class="container py-4">
    <div class="d-flex flex-wrap align-items-end justify-content-between gap-3 mb-3">
      <div>
        <h1 class="h3 mb-1">State Comparison Timeseries</h1>
        <div class="text-secondary">Multi-experiment dynamics across all Kira utilities</div>
      </div>
      <div class="metrics-pills">
        <span class="badge text-bg-dark" id="pill-exp">0 experiments</span>
        <span class="badge text-bg-primary" id="pill-metrics">0 metrics</span>
        <span class="badge text-bg-danger" id="pill-significant">0 significant deltas</span>
      </div>
    </div>

    <div class="card border-0 shadow-sm mb-4">
      <div class="card-body">
        <div class="row g-3 align-items-end">
          <div class="col-md-8">
            <label class="form-label">Metric</label>
            <select id="metricSelect" class="form-select"></select>
          </div>
          <div class="col-md-4">
            <label class="form-label">Quick filter</label>
            <input id="metricFilter" class="form-control" placeholder="e.g. redox, ros, nadh, qc">
          </div>
        </div>
      </div>
    </div>

    <div class="card border-0 shadow-sm mb-4">
      <div class="card-body">
        <div id="linePlot" class="plot"></div>
      </div>
    </div>

    <div class="card border-0 shadow-sm mb-4">
      <div class="card-header bg-white fw-semibold">Top significant changes</div>
      <div class="card-body p-0">
        <div class="table-responsive">
          <table class="table table-hover table-sm align-middle mb-0" id="sigTable">
            <thead>
              <tr>
                <th>Experiment</th>
                <th>Metric</th>
                <th class="text-end">Previous</th>
                <th class="text-end">Current</th>
                <th class="text-end">Delta</th>
              </tr>
            </thead>
            <tbody></tbody>
          </table>
        </div>
      </div>
    </div>
  </div>

  <script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>
  <script src="app.js"></script>
</body>
</html>
"#;

const STYLES_CSS: &str = r#"body {
  background:
    radial-gradient(1200px 300px at 10% -10%, rgba(13,110,253,0.18), transparent),
    radial-gradient(900px 260px at 90% 0%, rgba(220,53,69,0.16), transparent),
    #f6f8fb;
}

.plot {
  height: 430px;
}

.metrics-pills .badge {
  margin-left: 6px;
}
"#;

const APP_JS: &str = r#"async function main() {
  const r = await fetch('./data.json');
  const data = await r.json();

  const experiments = data.experiments || [];
  const metrics = data.metrics || {};
  const significant = data.significant || [];
  const metricKeys = Object.keys(metrics);

  document.getElementById('pill-exp').textContent = `${experiments.length} experiments`;
  document.getElementById('pill-metrics').textContent = `${metricKeys.length} metrics`;
  document.getElementById('pill-significant').textContent = `${significant.length} significant deltas`;

  const select = document.getElementById('metricSelect');
  const filter = document.getElementById('metricFilter');

  function refillSelect(query) {
    const q = (query || '').trim().toLowerCase();
    const current = select.value;
    select.innerHTML = '';
    for (const key of metricKeys) {
      const label = metrics[key].label || key;
      if (q && !label.toLowerCase().includes(q) && !key.toLowerCase().includes(q)) continue;
      const opt = document.createElement('option');
      opt.value = key;
      opt.textContent = label;
      select.appendChild(opt);
    }
    if (select.options.length === 0) return;
    if ([...select.options].some(o => o.value === current)) {
      select.value = current;
    }
    renderMetric(select.value);
  }

  function fmt(v) {
    if (v === null || Number.isNaN(v)) return '—';
    return Number(v).toFixed(4);
  }

  function renderMetric(key) {
    if (!key || !metrics[key]) return;
    const item = metrics[key];
    const y = item.values || [];
    Plotly.newPlot('linePlot', [{
      x: experiments,
      y,
      mode: 'lines+markers',
      marker: { size: 8, color: '#0d6efd' },
      line: { width: 2, color: '#0d6efd' },
      hovertemplate: '%{x}<br>%{y:.4f}<extra></extra>'
    }], {
      margin: { l: 60, r: 30, t: 40, b: 80 },
      title: { text: item.label || key, font: { size: 16 } },
      xaxis: { tickangle: -35, automargin: true },
      yaxis: { zeroline: true, gridcolor: '#e9ecef' },
      paper_bgcolor: 'white',
      plot_bgcolor: 'white'
    }, { responsive: true, displaylogo: false });
  }

  const tbody = document.querySelector('#sigTable tbody');
  const top = [...significant].slice(0, 200);
  for (const row of top) {
    const tr = document.createElement('tr');
    tr.innerHTML = `
      <td>${row.experiment}</td>
      <td title="${row.metric_key}">${row.metric_label}</td>
      <td class="text-end">${fmt(row.previous)}</td>
      <td class="text-end">${fmt(row.current)}</td>
      <td class="text-end fw-semibold ${Math.abs(row.delta) >= 0.25 ? 'text-danger' : 'text-warning'}">${fmt(row.delta)}</td>
    `;
    tbody.appendChild(tr);
  }

  select.addEventListener('change', () => renderMetric(select.value));
  filter.addEventListener('input', () => refillSelect(filter.value));

  refillSelect('');
}

main().catch(err => {
  console.error(err);
  document.body.insertAdjacentHTML('beforeend', `<pre class="m-3 p-3 bg-danger-subtle border border-danger-subtle rounded">Failed to load report data: ${String(err)}</pre>`);
});
"#;
