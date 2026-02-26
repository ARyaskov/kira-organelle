pub const CSS: &str = r#"
:root{
  --bg:#f6f7fb;
  --card:#ffffff;
  --text:#0f172a;
  --muted:#475569;
  --line:#d8deea;
  --adaptive:#1d4ed8;
  --transition:#d97706;
  --resistant:#dc2626;
}
*{box-sizing:border-box}
body{
  margin:0;
  font-family:ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
  color:var(--text);
  background:linear-gradient(180deg,#eef2ff 0%,var(--bg) 220px);
}
.wrap{max-width:1280px;margin:0 auto;padding:24px}
h1{margin:0 0 8px;font-size:30px}
h2{margin:0 0 14px;font-size:20px}
p,li{color:var(--muted)}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:16px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:12px}
.section{margin-top:18px}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}
.toolbar{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin-bottom:10px}
select,input[type=checkbox]{transform:translateY(1px)}
.chart{width:100%;height:300px;border:1px solid var(--line);border-radius:10px;background:#fff}
canvas.chart{display:block}
.legend{display:flex;gap:12px;flex-wrap:wrap}
.dot{display:inline-block;width:10px;height:10px;border-radius:999px;margin-right:6px}
.tooltip{
  position:fixed;
  pointer-events:none;
  background:#0b1222;
  color:#fff;
  padding:6px 8px;
  border-radius:8px;
  font-size:12px;
  opacity:0;
  transform:translate(8px,8px);
  white-space:nowrap;
}
.dimmed .chart{opacity:.35}
.decision-dot{cursor:pointer;stroke:#0b1222;stroke-width:1.2}
.decision-dot:focus{outline:none;stroke-width:2.5}
.decision-polyline{fill:none;stroke:#0b1222;stroke-width:1.5;stroke-dasharray:6 4}
.decision-legend{display:flex;gap:10px;flex-wrap:wrap;margin:8px 0}
"#;

pub const JS: &str = r#"
(function(){
const dataEl=document.getElementById('fii-data');
const report=JSON.parse(dataEl.textContent);
const decisionData=JSON.parse((document.getElementById('decision-data')||{textContent:'{}'}).textContent||'{}');
const decisionStabilityData=JSON.parse((document.getElementById('decision-stability-data')||{textContent:'{}'}).textContent||'{}');
const phaseData=JSON.parse((document.getElementById('phase-portrait-data')||{textContent:'{}'}).textContent||'{}');
const iliData=JSON.parse((document.getElementById('ili-data')||{textContent:'{}'}).textContent||'{}');
const caiData=JSON.parse((document.getElementById('cai-data')||{textContent:'{}'}).textContent||'{}');
const priData=JSON.parse((document.getElementById('pri-data')||{textContent:'{}'}).textContent||'{}');
const cocsData=JSON.parse((document.getElementById('cocs-data')||{textContent:'{}'}).textContent||'{}');
const dciData=JSON.parse((document.getElementById('dci-data')||{textContent:'{}'}).textContent||'{}');
let iliBySample=new Map();
const caiBySample=new Map((caiData.samples||[]).map(x=>[(x.sample_label||x.label||''),x]));
const priBySample=new Map((priData.samples||[]).map(x=>[(x.sample_label||x.label||''),x]));
const cocsBySample=new Map((cocsData.samples||[]).map(x=>[(x.sample_label||x.label||''),x]));
const dciBySample=new Map((dciData.samples||[]).map(x=>[(x.sample_label||x.label||''),x]));
const tip=document.getElementById('tooltip');
const pct=v=>(v*100).toFixed(1)+'%';
const fmt=v=>Number(v).toFixed(3);

function showTip(ev,text){
  const maxW=260;
  const maxH=80;
  const x=Math.min(window.innerWidth-maxW-8,Math.max(8,ev.clientX+10));
  const y=Math.min(window.innerHeight-maxH-8,Math.max(8,ev.clientY+10));
  tip.style.opacity='1';
  tip.style.left=x+'px';
  tip.style.top=y+'px';
  tip.textContent=text;
}
function hideTip(){tip.style.opacity='0'}
document.addEventListener('keydown',(ev)=>{if(ev.key==='Escape')hideTip();});

const cards=document.getElementById('overview-cards');
report.samples.forEach(s=>{
  const d=document.createElement('div');
  d.className='card';
  d.innerHTML=`<h3>${s.label}</h3>
  <div class="mono">${s.path}</div>
  <div>N cells: <b>${s.n_cells}</b></div>
  <div>FII mean/median: <b>${fmt(s.summary.mean)}</b> / <b>${fmt(s.summary.median)}</b></div>
  <div>Regimes: A <b>${pct(s.summary.adaptive)}</b>, T <b>${pct(s.summary.transition)}</b>, R <b>${pct(s.summary.resistant)}</b></div>
  <div>Low-confidence: <b>${pct(s.summary.low_confidence_fraction)}</b></div>`;
  cards.appendChild(d);
});

const modeSel=document.getElementById('dist-mode');
const distCanvas=document.getElementById('dist-canvas');
const regimeCanvas=document.getElementById('regime-canvas');
const sampleSel=document.getElementById('sample-select');
const regimeSel=document.getElementById('regime-filter');
const highConfOnly=document.getElementById('high-conf-only');
const sampleSelTraj=document.getElementById('traj-mode');
const overlayToggle=document.getElementById('show-decision-overlay');
const trajectoryToggle=document.getElementById('show-decision-trajectory');
const dimBaseToggle=document.getElementById('dim-base-plots');

report.samples.forEach((s,i)=>{
  const o=document.createElement('option');
  o.value=String(i); o.textContent=s.label;
  sampleSel.appendChild(o);
});

function setupCanvas(c){
  const dpr=Math.max(1,window.devicePixelRatio||1);
  const w=c.clientWidth||900,h=c.clientHeight||300;
  c.width=Math.floor(w*dpr); c.height=Math.floor(h*dpr);
  const ctx=c.getContext('2d');
  ctx.setTransform(dpr,0,0,dpr,0,0);
  return {ctx,w,h};
}

function drawDistribution(){
  const {ctx,w,h}=setupCanvas(distCanvas);
  ctx.clearRect(0,0,w,h);
  ctx.fillStyle='#fff';ctx.fillRect(0,0,w,h);
  ctx.strokeStyle='#d8deea';ctx.strokeRect(0,0,w,h);
  const mode=modeSel.value;
  const maxBin=Math.max(...report.samples.flatMap(s=>s.fii_histogram),1);
  if(mode==='overlay'){
    const colors=['#1d4ed8','#d97706','#dc2626','#0f766e','#9333ea','#be123c','#155e75'];
    report.samples.forEach((s,si)=>{
      ctx.beginPath();
      s.fii_histogram.forEach((v,i)=>{
        const x=(i/(s.fii_histogram.length-1))*w;
        const y=h-(v/maxBin)*(h-24)-12;
        if(i===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
      });
      ctx.strokeStyle=colors[si%colors.length];
      ctx.lineWidth=2;ctx.stroke();
    });
  }else{
    const rows=report.samples.length;
    const fh=Math.max(40,h/rows);
    report.samples.forEach((s,ri)=>{
      const y0=ri*fh;
      const localMax=Math.max(...s.fii_histogram,1);
      ctx.fillStyle='#0f172a';ctx.fillText(s.label,6,y0+12);
      s.fii_histogram.forEach((v,i)=>{
        const x=(i/s.fii_histogram.length)*w;
        const bw=w/s.fii_histogram.length;
        const bh=(v/localMax)*(fh-18);
        ctx.fillStyle='#94a3b8';
        ctx.fillRect(x,y0+fh-bh,bw-1,bh);
      });
    });
  }
}

function drawRegimes(){
  const {ctx,w,h}=setupCanvas(regimeCanvas);
  ctx.clearRect(0,0,w,h);
  const n=report.samples.length;
  const gap=10;
  const bw=(w-gap*(n+1))/Math.max(n,1);
  report.samples.forEach((s,i)=>{
    const x=gap+i*(bw+gap);
    let y=h-20;
    const segs=[
      ['adaptive',s.summary.adaptive,'#1d4ed8'],
      ['transition',s.summary.transition,'#d97706'],
      ['resistant',s.summary.resistant,'#dc2626']
    ];
    segs.forEach(([name,val,color])=>{
      const hh=(h-50)*val;
      y-=hh;
      ctx.fillStyle=color;ctx.fillRect(x,y,bw,hh);
    });
    ctx.fillStyle='#0f172a';
    ctx.save();ctx.translate(x+bw/2,h-6);ctx.rotate(-0.4);
    ctx.fillText(s.label,-24,0);ctx.restore();
  });
}

function heatmapFor(component){
  const s=report.samples[Number(sampleSel.value)||0];
  const entry=s.coupling.find(c=>c.component===component);
  const reg=regimeSel.value;
  const high=highConfOnly.checked;
  const key=(reg==='all'? 'all':reg)+(high?'_high_conf':'');
  return entry[key];
}

function drawCoupling(canvasId, component){
  const c=document.getElementById(canvasId);
  const {ctx,w,h}=setupCanvas(c);
  const payload=heatmapFor(component);
  const binsX=50,binsY=50;
  const max=Math.max(...payload.counts,1);
  for(let y=0;y<binsY;y++){
    for(let x=0;x<binsX;x++){
      const idx=y*binsX+x;
      const v=payload.counts[idx];
      const a=Math.sqrt(v/max);
      const px=(x/binsX)*w;
      const py=h-((y+1)/binsY)*h;
      ctx.fillStyle=`rgba(12,74,110,${a})`;
      ctx.fillRect(px,py,w/binsX+0.3,h/binsY+0.3);
    }
  }
  c.onmousemove=(ev)=>{
    const rect=c.getBoundingClientRect();
    const x=Math.floor(((ev.clientX-rect.left)/rect.width)*binsX);
    const y=binsY-1-Math.floor(((ev.clientY-rect.top)/rect.height)*binsY);
    if(x<0||x>=binsX||y<0||y>=binsY){hideTip();return;}
    const idx=y*binsX+x;
    const count=payload.counts[idx]||0;
    const mx=payload.mean_x[idx]||0;
    const my=payload.mean_y[idx]||0;
    showTip(ev,`${component}: n=${count}, mean x=${fmt(mx)}, mean FII=${fmt(my)}`);
  };
  c.onmouseleave=hideTip;
}

function drawTrajectory(){
  const c=document.getElementById('traj-canvas');
  const {ctx,w,h}=setupCanvas(c);
  const ordered=[...report.samples].sort((a,b)=>a.order-b.order);
  if(ordered.length<2){ctx.fillText('Need at least 2 ordered samples',12,20);return;}
  const mode=sampleSelTraj.value;
  const values=ordered.map(s=>mode==='mean'?s.summary.mean:s.summary.median);
  const values2=ordered.map(s=>s.summary.resistant);
  const drawLine=(vals,color,yScale)=>{ctx.beginPath();vals.forEach((v,i)=>{const x=(i/(vals.length-1))*(w-40)+20;const y=h-20-v*yScale;if(i===0)ctx.moveTo(x,y);else ctx.lineTo(x,y)});ctx.strokeStyle=color;ctx.lineWidth=2;ctx.stroke();};
  drawLine(values,'#1d4ed8',h-40);
  drawLine(values2,'#dc2626',h-40);
  ordered.forEach((s,i)=>{const x=(i/(ordered.length-1))*(w-40)+20;ctx.fillStyle='#0f172a';ctx.fillText(s.label,x-14,h-4)});
}

function tierColor(t){
  if(t==='STABLE_ADAPTIVE') return '#15803d';
  if(t==='TRANSITION_RISK') return '#d97706';
  if(t==='PRE_RESISTANT') return '#f97316';
  if(t==='FIXED_RESISTANT') return '#dc2626';
  return '#64748b';
}
function tierShape(t,x,y,size){
  if(t==='STABLE_ADAPTIVE') return `M ${x} ${y-size} L ${x+size} ${y} L ${x} ${y+size} L ${x-size} ${y} Z`;
  if(t==='TRANSITION_RISK') return `M ${x-size} ${y-size} H ${x+size} V ${y+size} H ${x-size} Z`;
  if(t==='PRE_RESISTANT') return `M ${x} ${y-size} L ${x+size} ${y+size} H ${x-size} Z`;
  return `M ${x-size} ${y-size} L ${x+size} ${y+size} M ${x+size} ${y-size} L ${x-size} ${y+size}`;
}
function clearSvg(el){while(el.firstChild)el.removeChild(el.firstChild);}
function addPoint(svg,x,y,sample,i,total){
  const g=document.createElementNS('http://www.w3.org/2000/svg','g');
  g.setAttribute('tabindex','0');
  g.setAttribute('role','button');
  g.setAttribute('aria-label',`${sample.sample_label||sample.label||'sample'} ${sample.decision_tier||sample.tier||''} confidence ${(sample.confidence||0).toFixed(2)}`);
  g.dataset.idx=String(i);
  const tier=sample.decision_tier||sample.tier||'UNKNOWN';
  const conf=Math.max(0.15,Math.min(1,sample.confidence||0));
  if(tier==='FIXED_RESISTANT'){
    const l1=document.createElementNS('http://www.w3.org/2000/svg','line');
    l1.setAttribute('x1',String(x-6));l1.setAttribute('y1',String(y-6));l1.setAttribute('x2',String(x+6));l1.setAttribute('y2',String(y+6));
    l1.setAttribute('stroke',tierColor(tier));l1.setAttribute('stroke-width','2.2');l1.setAttribute('opacity',String(conf));l1.classList.add('decision-dot');
    const l2=document.createElementNS('http://www.w3.org/2000/svg','line');
    l2.setAttribute('x1',String(x+6));l2.setAttribute('y1',String(y-6));l2.setAttribute('x2',String(x-6));l2.setAttribute('y2',String(y+6));
    l2.setAttribute('stroke',tierColor(tier));l2.setAttribute('stroke-width','2.2');l2.setAttribute('opacity',String(conf));l2.classList.add('decision-dot');
    g.appendChild(l1);g.appendChild(l2);
  }else{
    const p=document.createElementNS('http://www.w3.org/2000/svg','path');
    p.setAttribute('d',tierShape(tier,x,y,6));
    p.setAttribute('fill',tierColor(tier));
    p.setAttribute('opacity',String(conf));
    p.classList.add('decision-dot');
    g.appendChild(p);
  }
  const ili=iliBySample.get(sample.sample_label||sample.label||'');
  const iliTxt=ili?`${ili.ili||ili.ILI_category||'UNRESOLVED'} (${Number(ili.confidence||ili.ILI_confidence||0).toFixed(2)})`:'-';
  const cai=caiBySample.get(sample.sample_label||sample.label||'');
  const caiRaw=Number((sample.cai_context??sample.cai??(cai? (cai.CAI??cai.cai):NaN)));
  const caiTxt=Number.isFinite(caiRaw)?caiRaw.toFixed(2):'-';
  const pri=priBySample.get(sample.sample_label||sample.label||'');
  const priRaw=Number((sample.pri_context??sample.pri??(pri? (pri.PRI??pri.pri):NaN)));
  const priTxt=Number.isFinite(priRaw)?priRaw.toFixed(2):'-';
  const cocs=cocsBySample.get(sample.sample_label||sample.label||'');
  const cocsRaw=Number((sample.cocs_context??sample.cocs??(cocs? (cocs.COCS_global??cocs.cocs_global):NaN)));
  const cocsTxt=Number.isFinite(cocsRaw)?cocsRaw.toFixed(2):'-';
  const dci=dciBySample.get(sample.sample_label||sample.label||'');
  const dciRaw=Number((sample.dci_context??sample.dci??(dci? (dci.DCI??dci.dci):NaN)));
  const dciTxt=Number.isFinite(dciRaw)?dciRaw.toFixed(2):'-';
  const focusShow=(ev)=>showTip(ev,`Sample: ${sample.sample_label||sample.label} | Tier: ${tier} | Confidence: ${(sample.confidence||0).toFixed(2)} | PRI: ${priTxt} | CAI: ${caiTxt} | COCS: ${cocsTxt} | DCI: ${dciTxt} | ILI: ${iliTxt} | Drivers: ${(sample.drivers||[]).join(', ')||'-'}`);
  g.addEventListener('mousemove',focusShow);
  g.addEventListener('focus',focusShow);
  g.addEventListener('mouseleave',hideTip);
  g.addEventListener('blur',hideTip);
  g.addEventListener('keydown',(ev)=>{
    if(ev.key==='ArrowRight'){const n=svg.querySelector(`[data-idx="${Math.min(total-1,i+1)}"]`);if(n)n.focus();}
    if(ev.key==='ArrowLeft'){const p=svg.querySelector(`[data-idx="${Math.max(0,i-1)}"]`);if(p)p.focus();}
    if(ev.key==='Escape'){hideTip();}
  });
  svg.appendChild(g);
}
function drawDecisionOverlay(){
  const samples=(decisionData.samples||[]).slice().sort((a,b)=>(a.order_rank??a.rank??0)-(b.order_rank??b.rank??0));
  iliBySample=new Map((iliData.samples||[]).map(x=>[(x.sample_label||x.label||''),x]));
  const timelineSvg=document.getElementById('decision-timeline-svg');
  const phase1=document.getElementById('decision-phase-fii-vel');
  const phase2=document.getElementById('decision-phase-vel-acc');
  [timelineSvg,phase1,phase2].forEach(clearSvg);
  if(!overlayToggle.checked||samples.length===0) return;

  const W=1200,H=260,margin=40;
  const points=samples.map((s,i)=>({x:margin+(i*Math.max(1,(W-2*margin))/Math.max(1,samples.length-1)),y:H-40-(Math.max(0,Math.min(1,s.confidence||0))*(H-80)),s}));
  if(trajectoryToggle.checked&&points.length>1){
    const p=document.createElementNS('http://www.w3.org/2000/svg','polyline');
    p.setAttribute('points',points.map(v=>`${v.x},${v.y}`).join(' '));
    p.classList.add('decision-polyline');timelineSvg.appendChild(p);
  }
  points.forEach((v,i)=>{addPoint(timelineSvg,v.x,v.y,{...v.s,sample_label:v.s.sample_label||v.s.label},i,points.length);
    const tx=document.createElementNS('http://www.w3.org/2000/svg','text');tx.setAttribute('x',String(v.x-16));tx.setAttribute('y',String(H-12));tx.setAttribute('font-size','11');tx.textContent=v.s.sample_label||v.s.label||`S${i}`;timelineSvg.appendChild(tx);});

  const phaseSamples=(phaseData.samples||[]).slice().sort((a,b)=>(a.rank??0)-(b.rank??0));
  if(phaseSamples.length===0) return;
  const velVals=phaseSamples.map(s=>s.velocity).filter(v=>typeof v==='number');
  const accVals=phaseSamples.map(s=>s.acceleration).filter(v=>typeof v==='number');
  const velMin=velVals.length?Math.min(...velVals):-1,velMax=velVals.length?Math.max(...velVals):1;
  const accMin=accVals.length?Math.min(...accVals):-1,accMax=accVals.length?Math.max(...accVals):1;
  const pad=(a,b)=>{const s=Math.max(0.05,(b-a)*0.1);return [a-s,b+s];};
  const [v0,v1]=pad(velMin,velMax),[a0,a1]=pad(accMin,accMax);
  const scale=(v,min,max,size)=>30+((v-min)/(max-min||1))*(size-60);
  const H2=320,W2=580;
  const p1=phaseSamples.filter(s=>typeof s.velocity==='number').map((s,i)=>({x:scale(s.fii_mean??0,0,1,W2),y:H2-scale(s.velocity,v0,v1,H2),s,idx:i}));
  const p2=phaseSamples.filter(s=>typeof s.velocity==='number'&&typeof s.acceleration==='number').map((s,i)=>({x:scale(s.velocity,v0,v1,W2),y:H2-scale(s.acceleration,a0,a1,H2),s,idx:i}));
  if(trajectoryToggle.checked&&p1.length>1){const l=document.createElementNS('http://www.w3.org/2000/svg','polyline');l.setAttribute('points',p1.map(p=>`${p.x},${p.y}`).join(' '));l.classList.add('decision-polyline');phase1.appendChild(l);}
  if(trajectoryToggle.checked&&p2.length>1){const l=document.createElementNS('http://www.w3.org/2000/svg','polyline');l.setAttribute('points',p2.map(p=>`${p.x},${p.y}`).join(' '));l.classList.add('decision-polyline');phase2.appendChild(l);}
  p1.forEach((p,i)=>{const d=(decisionData.samples||[]).find(x=>(x.sample_label||x.label)===(p.s.label||p.s.sample_label))||{};addPoint(phase1,p.x,p.y,{...d,sample_label:p.s.label||p.s.sample_label,tier:d.decision_tier||d.tier,confidence:d.confidence||0,drivers:d.drivers||[]},i,p1.length);
    const ili=iliBySample.get(p.s.label||p.s.sample_label||'');
    const tier=d.decision_tier||d.tier||'';
    if(ili&&tier&&(tier==='PRE_RESISTANT'||tier==='FIXED_RESISTANT')){
      const t=document.createElementNS('http://www.w3.org/2000/svg','text');
      t.setAttribute('x',String(p.x+7));t.setAttribute('y',String(p.y-8));t.setAttribute('font-size','10');
      t.setAttribute('fill','#0f172a');t.textContent=(ili.ili||ili.ILI_category||'').replace('_DRIVEN','');
      phase1.appendChild(t);
    }});
  p2.forEach((p,i)=>{const d=(decisionData.samples||[]).find(x=>(x.sample_label||x.label)===(p.s.label||p.s.sample_label))||{};addPoint(phase2,p.x,p.y,{...d,sample_label:p.s.label||p.s.sample_label,tier:d.decision_tier||d.tier,confidence:d.confidence||0,drivers:d.drivers||[]},i,p2.length);
    const ili=iliBySample.get(p.s.label||p.s.sample_label||'');
    const tier=d.decision_tier||d.tier||'';
    if(ili&&tier&&(tier==='PRE_RESISTANT'||tier==='FIXED_RESISTANT')){
      const t=document.createElementNS('http://www.w3.org/2000/svg','text');
      t.setAttribute('x',String(p.x+7));t.setAttribute('y',String(p.y-8));t.setAttribute('font-size','10');
      t.setAttribute('fill','#0f172a');t.textContent=(ili.ili||ili.ILI_category||'').replace('_DRIVEN','');
      phase2.appendChild(t);
    }});
}

function applyDimBase(){
  const baseIds=['section-distribution','section-regime','section-coupling'];
  baseIds.forEach(id=>{const el=document.getElementById(id);if(!el)return;el.classList.toggle('dimmed',dimBaseToggle.checked);});
}

modeSel.onchange=drawDistribution;
sampleSel.onchange=()=>{drawCoupling('heatmap-mito','mitochondrial');drawCoupling('heatmap-translation','translation');drawCoupling('heatmap-splice','splice');};
regimeSel.onchange=()=>{drawCoupling('heatmap-mito','mitochondrial');drawCoupling('heatmap-translation','translation');drawCoupling('heatmap-splice','splice');};
highConfOnly.onchange=()=>{drawCoupling('heatmap-mito','mitochondrial');drawCoupling('heatmap-translation','translation');drawCoupling('heatmap-splice','splice');};
sampleSelTraj.onchange=drawTrajectory;
overlayToggle.onchange=drawDecisionOverlay;
trajectoryToggle.onchange=drawDecisionOverlay;
dimBaseToggle.onchange=applyDimBase;
window.addEventListener('resize',()=>{drawDistribution();drawRegimes();drawCoupling('heatmap-mito','mitochondrial');drawCoupling('heatmap-translation','translation');drawCoupling('heatmap-splice','splice');drawTrajectory();drawDecisionOverlay();});

drawDistribution();
drawRegimes();
drawCoupling('heatmap-mito','mitochondrial');
drawCoupling('heatmap-translation','translation');
drawCoupling('heatmap-splice','splice');
drawTrajectory();
drawDecisionOverlay();
applyDimBase();
})();
"#;
