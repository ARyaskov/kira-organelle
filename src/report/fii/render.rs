use serde::Serialize;
use serde_json::Value;

use super::assets::{CSS, JS};
use super::read::SampleReportData;

#[derive(Debug, Clone, Serialize)]
pub struct ReportPayload {
    pub title: String,
    pub generated_at: String,
    pub version: String,
    pub git_commit: Option<String>,
    pub samples: Vec<SampleReportData>,
    pub decision_data: Value,
    pub decision_stability_data: Value,
    pub phase_portrait_data: Value,
    pub ili_data: Value,
    pub cai_data: Value,
    pub pri_data: Value,
    pub cocs_data: Value,
    pub dci_data: Value,
}

pub fn render_html(payload: &ReportPayload) -> Result<String, String> {
    let json = serde_json::to_string(payload).map_err(|e| format!("PARSE payload JSON: {e}"))?;
    let commit = payload.git_commit.as_deref().unwrap_or("unknown");
    let decision_json = serde_json::to_string(&payload.decision_data)
        .map_err(|e| format!("PARSE decision-data JSON: {e}"))?;
    let decision_stability_json = serde_json::to_string(&payload.decision_stability_data)
        .map_err(|e| format!("PARSE decision-stability-data JSON: {e}"))?;
    let phase_json = serde_json::to_string(&payload.phase_portrait_data)
        .map_err(|e| format!("PARSE phase-portrait-data JSON: {e}"))?;
    let ili_json = serde_json::to_string(&payload.ili_data)
        .map_err(|e| format!("PARSE ili-data JSON: {e}"))?;
    let cai_json = serde_json::to_string(&payload.cai_data)
        .map_err(|e| format!("PARSE cai-data JSON: {e}"))?;
    let pri_json = serde_json::to_string(&payload.pri_data)
        .map_err(|e| format!("PARSE pri-data JSON: {e}"))?;
    let cocs_json = serde_json::to_string(&payload.cocs_data)
        .map_err(|e| format!("PARSE cocs-data JSON: {e}"))?;
    let dci_json = serde_json::to_string(&payload.dci_data)
        .map_err(|e| format!("PARSE dci-data JSON: {e}"))?;

    Ok(format!(
        r#"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{title}</title>
<style>{css}</style>
</head>
<body>
<div class="wrap">
  <div class="card">
    <h1>{title}</h1>
    <p>Generated at: <span class="mono">{generated}</span></p>
    <p>kira-organelle version: <span class="mono">{version}</span> | git commit: <span class="mono">{commit}</span></p>
    <div id="sample-list" class="mono"></div>
  </div>

  <div class="section">
    <h2>Overview Cards</h2>
    <div id="overview-cards" class="grid"></div>
  </div>

  <div class="section" id="section-distribution">
    <h2>Distribution View</h2>
    <div class="toolbar">
      <label>Mode
        <select id="dist-mode">
          <option value="overlay">Overlay</option>
          <option value="facet">Facet by sample</option>
        </select>
      </label>
      <span class="mono">Stable bins: 50 over [0,1]</span>
    </div>
    <canvas id="dist-canvas" class="chart" width="1200" height="300"></canvas>
  </div>

  <div class="section" id="section-regime">
    <h2>Regime Composition</h2>
    <div class="legend">
      <span><i class="dot" style="background:var(--adaptive)"></i>Adaptive</span>
      <span><i class="dot" style="background:var(--transition)"></i>Transition</span>
      <span><i class="dot" style="background:var(--resistant)"></i>Resistant</span>
    </div>
    <canvas id="regime-canvas" class="chart" width="1200" height="300"></canvas>
  </div>

  <div class="section" id="section-coupling">
    <h2>Component Coupling View</h2>
    <div class="toolbar">
      <label>Sample
        <select id="sample-select"></select>
      </label>
      <label>Regime
        <select id="regime-filter">
          <option value="all">All</option>
          <option value="adaptive">Adaptive</option>
          <option value="transition">Transition</option>
          <option value="resistant">Resistant</option>
        </select>
      </label>
      <label><input type="checkbox" id="high-conf-only"/> LOW_CONFIDENCE = false only</label>
    </div>
    <div class="grid">
      <div class="card"><h3>FII vs mito score</h3><canvas id="heatmap-mito" class="chart" width="380" height="280"></canvas></div>
      <div class="card"><h3>FII vs translation score</h3><canvas id="heatmap-translation" class="chart" width="380" height="280"></canvas></div>
      <div class="card"><h3>FII vs splice score</h3><canvas id="heatmap-splice" class="chart" width="380" height="280"></canvas></div>
    </div>
  </div>

  <section class="section" id="decision-overlay">
    <h2>Decision Overlay</h2>
    <div class="toolbar">
      <label><input type="checkbox" id="show-decision-overlay" checked/> Show decision overlay</label>
      <label><input type="checkbox" id="show-decision-trajectory" checked/> Show decision trajectory</label>
      <label><input type="checkbox" id="dim-base-plots"/> Dim base plots</label>
    </div>
    <div class="grid">
      <div class="card">
        <h3>Decision Timeline Overlay</h3>
        <svg id="decision-timeline-svg" class="chart" viewBox="0 0 1200 260" role="img" aria-label="Decision timeline overlay"></svg>
      </div>
      <div class="card">
        <h3>Phase Projection: FII vs Velocity</h3>
        <svg id="decision-phase-fii-vel" class="chart" viewBox="0 0 580 320" role="img" aria-label="Decision projection FII versus Velocity"></svg>
      </div>
      <div class="card">
        <h3>Phase Projection: Velocity vs Acceleration</h3>
        <svg id="decision-phase-vel-acc" class="chart" viewBox="0 0 580 320" role="img" aria-label="Decision projection Velocity versus Acceleration"></svg>
      </div>
    </div>
  </section>

  <div class="section">
    <h2>Sample Trajectory Summary</h2>
    <div class="toolbar">
      <label>Line
        <select id="traj-mode">
          <option value="mean">Mean FII</option>
          <option value="median">Median FII</option>
        </select>
      </label>
      <span class="mono">Red line: Resistant fraction</span>
    </div>
    <canvas id="traj-canvas" class="chart" width="1200" height="260"></canvas>
  </div>
</div>

<script id="fii-data" type="application/json">{json}</script>
<script id="decision-data" type="application/json">{decision_json}</script>
<script id="decision-stability-data" type="application/json">{decision_stability_json}</script>
<script id="phase-portrait-data" type="application/json">{phase_json}</script>
<script id="ili-data" type="application/json">{ili_json}</script>
<script id="cai-data" type="application/json">{cai_json}</script>
<script id="pri-data" type="application/json">{pri_json}</script>
<script id="cocs-data" type="application/json">{cocs_json}</script>
<script id="dci-data" type="application/json">{dci_json}</script>
<div id="tooltip" class="tooltip"></div>
<script>{js}</script>
<script>
document.getElementById('sample-list').textContent =
  'Samples: ' + JSON.parse(document.getElementById('fii-data').textContent).samples.map(s => s.label).join(', ');
</script>
</body>
</html>
"#,
        title = escape_html(&payload.title),
        generated = escape_html(&payload.generated_at),
        version = escape_html(&payload.version),
        commit = escape_html(commit),
        css = CSS,
        js = JS,
        json = json,
        decision_json = decision_json,
        decision_stability_json = decision_stability_json,
        phase_json = phase_json,
        ili_json = ili_json,
        cai_json = cai_json,
        pri_json = pri_json,
        cocs_json = cocs_json,
        dci_json = dci_json
    ))
}

fn escape_html(raw: &str) -> String {
    raw.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}
