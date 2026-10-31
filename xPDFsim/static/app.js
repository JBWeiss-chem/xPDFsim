import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const $ = (id) => document.getElementById(id);
// opened structures; every entry keeps its CIF text, structure and calculated histograms
const S = {
  entries: [], active: null, conv: 'G', rad: 'xray', sigma: 0.1, showHist: false, hl: null, hlWidth: 0.1, hlMin: 0.1, hlRod: 0.2, hlDim: 0.38, dragMode: 'zoom', seen: new Set(),
  hiddenBy: { xray: new Set(), neutron: new Set() },
  showRef: false, xr: null, sc: [1, 1, 1],
};
// curves hidden in the plot, kept separately for X-rays and neutrons
Object.defineProperty(S, 'hidden', { get: () => S.hiddenBy[S.rad] });
let nextId = 1;
const PALETTE = ['#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f', '#edc948', '#b07aa1',
  '#ff9da7', '#9c755f', '#bab0ac', '#86bcb6', '#d37295', '#fabfd2', '#8cd17d', '#499894', '#f1ce63'];
const INK = '#1d2330';
const OUT_STEP = 0.02;          // Å, grid of the broadened curves (as in the command line tool)
const MIN_BAR_PX = 4;           // minimum width of a histogram bar on screen

// Real-space functions of Peterson, Olds, McDonnell & Page, J. Appl. Cryst. 54, 317 (2021), Table 4,
// with the corrigendum by Peterson & Keen, J. Appl. Cryst. 54, 1542 (2021). All are derived from the
// weighted radial distribution function R(r); b = Z·rₑ (X-ray scattering length at Q = 0) and
// ⟨b⟩² is in barn. Every function except N(r) is linear in R(r): X(r) = a(r)·R(r) + c(r)·W, where
// W = 1 for the total and the Faber–Ziman weight of the pair for a partial (Appendices A1 and F).
// 'lim' gives the low-r and high-r limits of Table 4.
const PI4 = 4 * Math.PI;
const CONV = {
  R: {
    name: 'R(r)', file: 'R_r', unit: 'atom Å⁻¹', legend: 'left',
    terms: () => [1, 0], lim: (r, rho) => [0, PI4 * rho * r * r],
    formula: 'R(r) = 1/N Σ<sub>μν</sub> b<sub>μ</sub>b<sub>ν</sub>/⟨b⟩² δ(r − r<sub>μν</sub>)', eq: 'eq. 15',
    note: 'weighted radial distribution function',
  },
  rho: {
    name: 'ρ(r)', file: 'rho_r', unit: 'atom Å⁻³', legend: 'right',
    terms: (r) => [1 / (PI4 * r * r), 0], lim: (r, rho) => [0, rho],
    formula: 'ρ(r) = R(r) / 4πr²', eq: 'eq. 16', note: 'density function',
  },
  g: {
    name: 'g(r)', file: 'g_r', unit: '', legend: 'right',
    terms: (r, rho) => [1 / (PI4 * rho * r * r), 0], lim: () => [0, 1],
    formula: 'g(r) = ρ(r)/ρ₀ = R(r) / 4πρ₀r²', eq: 'eq. 22', note: 'pair distribution function (recommended)',
  },
  GK: {
    name: 'G_K(r)', file: 'GK_r', unit: 'barn', legend: 'right',
    terms: (r, rho, b2) => [b2 / (PI4 * rho * r * r), -b2], lim: (r, rho, b2) => [-b2, 0],
    formula: 'G<sub>K</sub>(r) = ⟨b⟩² [g(r) − 1]', eq: 'eq. 26', note: 'total radial distribution function (Keen 2001)',
  },
  G: {
    name: 'G(r)', file: 'G_r', unit: 'atom Å⁻²', legend: 'right',
    terms: (r, rho) => [1 / r, -PI4 * rho * r], lim: (r, rho) => [-PI4 * rho * r, 0],
    formula: 'G(r) = 4πr [ρ(r) − ρ₀] = R(r)/r − 4πρ₀r', eq: 'eq. 19', note: 'reduced pair distribution function (recommended)',
  },
  D: {
    name: 'D(r)', file: 'D_r', unit: 'barn atom Å⁻²', legend: 'right',
    terms: (r, rho, b2) => [b2 / r, -PI4 * rho * b2 * r], lim: (r, rho, b2) => [-PI4 * rho * b2 * r, 0],
    formula: 'D(r) = ⟨b⟩² G(r)', eq: 'eq. 28', note: 'differential correlation function',
  },
  T: {
    name: 'T(r)', file: 'T_r', unit: 'barn atom Å⁻²', legend: 'left',
    terms: (r, rho, b2) => [b2 / r, 0], lim: (r, rho, b2) => [0, PI4 * rho * b2 * r],
    formula: 'T(r) = ⟨b⟩² R(r)/r = ⟨b⟩² [G(r) + 4πρ₀r]', eq: 'eqs. 29, 32', note: 'total correlation function',
  },
  N: {
    name: 'N(r)', file: 'N_r', unit: 'atom', legend: 'left', totalOnly: true,
    lim: (r, rho) => [0, (PI4 / 3) * rho * r ** 3],
    formula: 'N(r) = ∫₀<sup>r</sup> RDF(r′) dr′', eq: 'eq. 31', note: 'coordination number (unweighted, total only)',
  },
};
// names with a subscript (G_K) for HTML and SVG
const nameHtml = (n) => n.replace(/_(\w)/, '<sub>$1</sub>');
const nameSvg = (n) => n.replace(/_(\w)/, '<tspan baseline-shift="sub" font-size="75%">$1</tspan>');

// ------------------------------------------------------------------ API
async function api(path, body) {
  const r = await fetch(path, {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
  });
  const j = await r.json();
  if (j.error) throw new Error(j.error);
  return j;
}
function fail(e) { console.error(e); alert(e.message || e); }
// every text that comes from a structure file or the server is escaped before it is put into HTML
const ESC = { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;', '`': '&#96;' };
const esc = (v) => String(v).replace(/[&<>"'`]/g, (c) => ESC[c]);

// ------------------------------------------------------------------ loading
function calcSettings() { return { r_max: +$('rMax').value, dr: +$('binWidth').value }; }

async function openFiles(files) {
  for (const f of files) {
    try { await loadText(f.name, await f.text()); } catch (e) { fail(e); }
  }
}

// name: shown for files with a single structure (default: the file name without extension)
async function loadText(file, text, name = file.replace(/\.cif$/i, '')) {
  const placeholder = { id: nextId++, file, name, loading: true };
  S.entries.push(placeholder);
  renderEntries();
  try {
    const r = await api('/api/load', { text, ...calcSettings() });
    const added = r.structures.map((s) => ({
      id: nextId++, file, text, block: s.block, name: r.structures.length > 1 ? s.name : name,
      st: s.structure, pdf: s.pdf, warnings: r.warnings,
    }));
    S.entries.splice(S.entries.indexOf(placeholder), 1, ...added);
    activate(added[0]);
    return added[0];
  } catch (e) {
    S.entries.splice(S.entries.indexOf(placeholder), 1);
    renderEntries();
    throw e;
  }
}

function activate(entry) {
  S.active = entry;
  S.xr = null;
  S.hl = null;
  $('readout').hidden = true;
  if (entry) {
    $('rMax').value = entry.pdf.r_max;
    $('binWidth').value = entry.pdf.dr;
    // every structure keeps its own supercell; a new one shows one unit cell
    entry.sc = entry.sc || [1, 1, 1];
    S.sc = entry.sc;
  }
  showSupercell();
  // element pairs with hydrogen are hidden at first for X-rays, where H scatters weakly; for
  // neutrons, where H scatters strongly, they are shown
  if (entry) Object.keys(entry.pdf.partials).forEach((l) => {
    if (S.seen.has(l)) return;
    S.seen.add(l);
    if (l.split('-').some((el) => el === 'H' || el === 'D')) S.hiddenBy.xray.add(l);
  });
  renderEntries();
  renderStructInfo();
  renderCalcInfo();
  buildScene();
  frameCamera();
  renderCurveList();
  updatePdf(false);
}

function closeEntry(entry) {
  const k = S.entries.indexOf(entry);
  S.entries.splice(k, 1);
  if (S.active === entry) activate(S.entries[Math.min(k, S.entries.length - 1)] || null);
  else renderEntries();
}

function renderEntries() {
  const list = $('structList');
  list.hidden = !S.entries.length;
  list.replaceChildren(...S.entries.map((e) => {
    const item = document.createElement('div');
    item.className = `sitem${e === S.active ? ' active' : ''}${e.loading ? ' loading' : ''}`;
    item.innerHTML = `<button class="pick" title="${esc(e.file)}"><span class="mark${e.loading || e.calculating ? ' running' : ''}"></span>
      <span class="nm">${esc(e.name)}</span><span class="file">${esc(e.file)}</span></button>
      <button class="x" title="close">×</button>`;
    if (!e.loading) {
      item.querySelector('.pick').onclick = () => { if (e !== S.active) activate(e); };
      item.querySelector('.x').onclick = () => closeEntry(e);
    }
    return item;
  }));
}

function renderStructInfo() {
  const e = S.active;
  $('warnings').innerHTML = e ? e.warnings.map((w) => `<div class="warn">${esc(w)}</div>`).join('') : '';
  if (!e) { $('structInfo').innerHTML = 'Load a CIF file to start.'; $('structInfo').className = 'muted'; return; }
  const st = e.st, f = (v, n) => v.toFixed(n);
  const formula = esc(st.formula).replace(/(\d+(\.\d+)?)/g, '<sub>$1</sub>');
  $('structInfo').className = '';
  $('structInfo').innerHTML = `<div class="kv">
    <span>formula</span><span>${formula}</span>
    <span>space group</span><span>${st.space_group ? esc(st.space_group) : '–'}</span>
    <span>a, b, c</span><span class="num">${st.abc.map((v) => f(v, 3)).join(', ')} Å</span>
    <span>α, β, γ</span><span class="num">${st.angles.map((v) => f(v, 2)).join(', ')}°</span>
    <span>volume</span><span class="num">${f(st.volume, 1)} Å³</span>
    <span>sites</span><span class="num">${st.atoms.length} in the cell${st.ordered ? '' : ', <b>disordered</b>'}</span>
    ${e.source ? `<span>source</span><span>COD ${esc(e.source.cod)}; ${esc(e.source.reference)}
      (<a href="https://doi.org/${encodeURI(e.source.doi)}" target="_blank" rel="noopener">doi</a>)</span>` : ''}
  </div>`;
}

// ------------------------------------------------------------------ calculation
function calcChanged() {
  const e = S.active;
  if (!e) return false;
  const { r_max, dr } = calcSettings();
  return r_max !== e.pdf.r_max || dr !== e.pdf.dr;
}
function renderCalcInfo() {
  const e = S.active;
  $('calc').disabled = !e || e.calculating || !calcChanged();
  $('calc').textContent = e && calcChanged() ? 'Recalculate PDF' : 'Calculate PDF';
  $('progress').hidden = !(e && e.calculating);
  if (!e) { $('calcInfo').textContent = ''; return; }
  const p = e.pdf;
  $('calcInfo').innerHTML = `Calculated up to ${p.r_max} Å with ${p.dr} Å bins in ${p.seconds} s on the unit cell
    (${+p.n_atoms.toFixed(3)} atoms); ρ₀ = ${p.rho.toFixed(5)} Å⁻³.${calcChanged() ? '<br><b>The settings above differ: recalculate to apply them.</b>' : ''}`;
}
['rMax', 'binWidth'].forEach((id) => { $(id).oninput = renderCalcInfo; });

$('calc').onclick = async () => {
  const e = S.active;
  if (!e) return;
  e.calculating = true;
  renderCalcInfo(); renderEntries();
  try {
    const { pdf } = await api('/api/pdf', { text: e.text, block: e.block, ...calcSettings() });
    e.pdf = pdf;
    e.dist = e.Rtot = e.pairCache = null;
    if (e === S.active) {
      if (S.xr && S.xr[1] > pdf.r_max) S.xr = null;
      if (S.hl && S.hl.r > pdf.r_max) S.hl = null;
      updatePdf(false);
      updateHighlight();
    }
  } catch (err) { fail(err); } finally {
    e.calculating = false;
    renderEntries();
    if (e === S.active) renderCalcInfo();
  }
};

// ------------------------------------------------------------------ PDF maths
// histograms -> chosen convention -> broadened curves on the display grid (and the histograms themselves)
// weights of the chosen radiation. The histograms are calculated with X-ray weights; a neutron
// partial is the X-ray one times the ratio of the bb factors (b_i b_j / <b>^2) of the pair.
function radiationWeights(pdf, rad) {
  if (rad !== 'neutron') return { weights: pdf.weights, mean_b2: pdf.mean_b2, scale: () => 1 };
  const n = pdf.neutron;
  return { weights: n.weights, mean_b2: n.mean_b2, scale: (l) => n.weights[l].bb / pdf.weights[l].bb };
}
function transform(pdf, conv, rad) {
  const n = pdf.n, rho = pdf.rho;
  const r = Float64Array.from({ length: n }, (_, k) => pdf.r0 + k * pdf.dr);
  const labels = Object.keys(pdf.partials);
  if (CONV[conv].totalOnly) {
    // N(r): integral of the unweighted RDF(r) (eq. 14) up to the bin centres
    const rdf = new Float64Array(n);
    labels.forEach((l) => pdf.partials[l].forEach((y, k) => { rdf[k] += y / pdf.weights[l].bb; }));
    let sum = 0;
    const N = rdf.map((y) => { sum += y; return (sum - y / 2) * pdf.dr; });
    return { r, curves: [{ label: 'total', color: INK, total: true, y: N }] };
  }
  const { weights, mean_b2, scale } = radiationWeights(pdf, rad);
  const partials = Object.fromEntries(labels.map((l) => [l, Float64Array.from(pdf.partials[l], (y) => y * scale(l))]));
  const total = new Float64Array(n);
  labels.forEach((l) => partials[l].forEach((y, k) => { total[k] += y; }));
  const terms = Array.from(r, (x) => CONV[conv].terms(x, rho, mean_b2));
  const apply = (y, W) => Float64Array.from(y, (v, k) => terms[k][0] * v + terms[k][1] * W);
  const curves = labels.map((l, i) => ({ label: l, color: PALETTE[i % PALETTE.length], y: apply(partials[l], weights[l].W) }));
  curves.push({ label: 'total', color: INK, total: true, y: apply(total, 1) });
  return { r, curves };
}

function broaden(r0, dr, y, xs, sigma) {
  // sum of Gaussians centred at every bin (truncated at 6σ, where they are < 1e-7 of their peak)
  const out = new Float64Array(xs.length);
  const norm = dr / (sigma * Math.sqrt(2 * Math.PI)), inv = 1 / (2 * sigma * sigma), w = 6 * sigma;
  for (let i = 0; i < xs.length; i++) {
    const x = xs[i];
    const k0 = Math.max(0, Math.ceil((x - w - r0) / dr)), k1 = Math.min(y.length - 1, Math.floor((x + w - r0) / dr));
    let s = 0;
    for (let k = k0; k <= k1; k++) { const d = x - (r0 + k * dr); s += y[k] * Math.exp(-d * d * inv); }
    out[i] = s * norm;
  }
  return out;
}

function computeCurves(entry, conv, rad) {
  const pdf = entry.pdf;
  const { r, curves } = transform(pdf, conv, rad);
  // both layers are always calculated, so the 'show' choice only changes what is drawn
  const n = r.findIndex((v) => v > pdf.r_max);
  const hx = r.slice(0, n < 0 ? r.length : n);
  const hist = curves.map((c) => ({ ...c, y: c.y.slice(0, hx.length) }));
  const x = Float64Array.from({ length: Math.round(pdf.r_max / OUT_STEP) + 1 }, (_, i) => i * OUT_STEP);
  const broad = curves.map((c) => ({ ...c, y: broaden(pdf.r0, pdf.dr, c.y, x, S.sigma) }));
  // low-r and high-r limits (Table 4), drawn as dashed asymptotes
  const mean_b2 = radiationWeights(pdf, rad).mean_b2;
  const lims = Array.from(x, (v) => CONV[conv].lim(v, pdf.rho, mean_b2));
  const refs = [0, 1].map((k) => Float64Array.from(lims, (l) => l[k]));
  return { conv, rad, x, curves: broad, hx, hist, dr: pdf.dr, refs, rMax: pdf.r_max };
}

// ------------------------------------------------------------------ chart
let shown = null;          // data on screen (possibly an interpolated tween frame)
let current = null;        // target data of the current settings
let tween = null;
let plot = null;           // geometry of the last drawn chart (for mouse interaction)

function updatePdf(animate) {
  const e = S.active;
  // neutron functions need tabulated scattering lengths for all elements
  const nError = e && e.pdf.neutron && e.pdf.neutron.error;
  document.querySelectorAll('[data-rad="neutron"]').forEach((b) => { b.disabled = !!nError; b.title = nError || ''; });
  if (nError && S.rad === 'neutron') S.rad = 'xray';
  const cv = CONV[S.conv];
  const radNote = S.rad === 'neutron'
    ? `neutron: b = coherent scattering length${e ? `, ⟨b⟩² = ${e.pdf.neutron.mean_b2.toPrecision(4)} barn` : ''}`
    : `X-ray: b = Z·rₑ (Q = 0)${e ? `, ⟨b⟩² = ${e.pdf.mean_b2.toPrecision(4)} barn` : ''}`;
  $('formula').innerHTML = `${cv.formula}<br><span class="note">${cv.note}; ${cv.eq}${cv.unit ? `; unit: ${cv.unit}` : ''}<br>${radNote}</span>`;
  document.querySelectorAll('[data-conv]').forEach((b) => b.classList.toggle('on', b.dataset.conv === S.conv && b.dataset.rad === S.rad));
  $('pdfTitle').innerHTML = e ? `${S.rad === 'neutron' ? 'Neutron' : 'X-ray'} ${nameHtml(cv.name)} of ${esc(e.name)}` : 'Pair distribution function';
  renderCurveList();
  if (!e) {
    current = shown = null; tween = null;
    $('chart').innerHTML = '<p class="muted empty">The PDF of the loaded structure appears here.</p>';
    updateControlStates();
    return;
  }
  const next = computeCurves(e, S.conv, S.rad);
  if (animate && shown && shown.x.length === next.x.length && shown.hx.length === next.hx.length) {
    tween = { from: shown, to: next, t0: performance.now() };
  } else {
    tween = null;
    shown = next;
    drawChart(shown);
  }
  current = next;
  syncRange();
}

// smooth transition between conventions: curves morph, appearing/vanishing curves fade
const TWEEN_MS = 380;
function tweenFrame(now) {
  if (tween) {
    const t = Math.min(1, (now - tween.t0) / TWEEN_MS);
    const s = t < 0.5 ? 2 * t * t : 1 - (-2 * t + 2) ** 2 / 2;
    const { from, to } = tween;
    const lerp = (a, b) => a.map((v, i) => v + (b[i] - v) * s);
    const morph = (fromCurves, toCurves) => {
      const fromBy = new Map(fromCurves.map((c) => [c.label, c]));
      const toLabels = new Set(toCurves.map((c) => c.label));
      const out = toCurves.map((c) => {
        const f = fromBy.get(c.label);
        return f ? { ...c, y: lerp(f.y, c.y), alpha: 1 } : { ...c, alpha: s };
      });
      fromCurves.forEach((c) => { if (!toLabels.has(c.label)) out.push({ ...c, alpha: 1 - s }); });
      return out;
    };
    const yr0 = yRange(from), yr1 = yRange(to);
    shown = { ...to, curves: morph(from.curves, to.curves), hist: morph(from.hist, to.hist), refs: to.refs.map((y, k) => lerp(from.refs[k], y)), yr: [yr0[0] + (yr1[0] - yr0[0]) * s, yr0[1] + (yr1[1] - yr0[1]) * s] };
    drawChart(shown);
    if (t >= 1) { tween = null; shown = to; drawChart(shown); }
  }
  requestAnimationFrame(tweenFrame);
}
requestAnimationFrame(tweenFrame);

function xRange(data) { return S.xr || [0, data.rMax]; }
function visible(c) { return !S.hidden.has(c.label); }
// the layer on screen: the broadened function or the histogram (never both)
function mainLayer(data) {
  return S.showHist ? { x: data.hx, curves: data.hist } : { x: data.x, curves: data.curves };
}
function yRange(data) {
  const [x0, x1] = xRange(data);
  let lo = Infinity, hi = -Infinity;
  const scan = (xs, y) => xs.forEach((x, i) => { if (x >= x0 && x <= x1) { lo = Math.min(lo, y[i]); hi = Math.max(hi, y[i]); } });
  const layer = mainLayer(data);
  layer.curves.filter(visible).forEach((c) => scan(layer.x, c.y));
  if (S.showRef) data.refs.forEach((y) => scan(data.x, y));
  if (!isFinite(lo)) { lo = 0; hi = 1; }
  lo = Math.min(lo, 0); hi = Math.max(hi, 0);
  if (hi - lo < 1e-9) { hi += 1; }
  const pad = 0.05 * (hi - lo);
  return [lo - (lo < 0 ? pad : 0), hi + pad];
}

function niceTicks(lo, hi, n) {
  const raw = (hi - lo) / n, mag = 10 ** Math.floor(Math.log10(raw)), f = raw / mag;
  const step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag;
  const ticks = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-6; v += step) ticks.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  const dec = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return { ticks, fmt: (v) => v.toFixed(dec) };
}

function drawChart(data) {
  const box = $('chart');
  const W = box.clientWidth, H = box.clientHeight;
  if (W < 50 || H < 50) return;
  const ml = 64, mr = 16, mt = 14, mb = 42;
  const [x0, x1] = xRange(data), [y0, y1] = data.yr || yRange(data);
  const X = (v) => ml + ((v - x0) / (x1 - x0)) * (W - ml - mr);
  const Y = (v) => mt + ((y1 - v) / (y1 - y0)) * (H - mt - mb);
  plot = { ml, mr, mt, mb, W, H, x0, x1, X, Y };
  const font = 'font-family="-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, sans-serif"';
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" ${font} font-size="12">`;
  svg += `<rect width="${W}" height="${H}" fill="#ffffff"/>`;
  svg += `<defs><clipPath id="clip"><rect x="${ml}" y="${mt}" width="${W - ml - mr}" height="${H - mt - mb}"/></clipPath></defs>`;
  // grid and tick labels
  const xt = niceTicks(x0, x1, Math.max(3, Math.floor((W - ml - mr) / 80)));
  const yt = niceTicks(y0, y1, Math.max(3, Math.floor((H - mt - mb) / 45)));
  xt.ticks.forEach((v) => {
    svg += `<line x1="${X(v)}" x2="${X(v)}" y1="${mt}" y2="${H - mb}" stroke="#eef0f3"/>`;
    svg += `<text x="${X(v)}" y="${H - mb + 16}" text-anchor="middle" fill="#6b7385">${xt.fmt(v)}</text>`;
  });
  yt.ticks.forEach((v) => {
    svg += `<line x1="${ml}" x2="${W - mr}" y1="${Y(v)}" y2="${Y(v)}" stroke="#eef0f3"/>`;
    svg += `<text x="${ml - 6}" y="${Y(v) + 4}" text-anchor="end" fill="#6b7385">${yt.fmt(v)}</text>`;
  });
  if (y0 < 0 && y1 > 0) svg += `<line x1="${ml}" x2="${W - mr}" y1="${Y(0)}" y2="${Y(0)}" stroke="#b9c0cc"/>`;
  svg += `<rect x="${ml}" y="${mt}" width="${W - ml - mr}" height="${H - mt - mb}" fill="none" stroke="#c9ced8"/>`;
  // axis labels
  const conv = CONV[data.conv];
  svg += `<text x="${(ml + W - mr) / 2}" y="${H - 8}" text-anchor="middle" fill="${INK}" font-size="13"><tspan font-style="italic">r</tspan> / Å</text>`;
  const yl = nameSvg(conv.unit ? `${conv.name} / ${conv.unit}` : conv.name);
  svg += `<text transform="translate(16 ${(mt + H - mb) / 2}) rotate(-90)" text-anchor="middle" fill="${INK}" font-size="13">${yl}</text>`;
  // curves: only the points inside the visible range (plus one on each side)
  const visibleRange = (xs) => {
    const i0 = Math.max(0, xs.findIndex((v) => v >= x0) - 1);
    const i1 = xs.findIndex((v) => v > x1);
    return [i0, i1 < 0 ? xs.length - 1 : i1];
  };
  const path = (xs, y) => {
    const [i0, i1] = visibleRange(xs);
    let d = '';
    for (let i = i0; i <= i1; i++) d += `${i === i0 ? 'M' : 'L'}${X(xs[i]).toFixed(1)} ${Y(y[i]).toFixed(1)}`;
    return d;
  };
  // histogram: one filled bar per bin, at least MIN_BAR_PX wide so that narrow bins stay visible
  const bars = (xs, y) => {
    const [i0, i1] = visibleRange(xs), y0 = Y(0);
    const half = Math.max(X(data.dr) - X(0), MIN_BAR_PX) / 2;
    let d = '';
    for (let i = i0; i <= i1; i++) {
      const top = Y(y[i]);
      if (Math.abs(top - y0) < 0.25) continue;
      const cx = X(xs[i]);
      d += `M${(cx - half).toFixed(1)} ${y0.toFixed(1)}V${top.toFixed(1)}H${(cx + half).toFixed(1)}V${y0.toFixed(1)}Z`;
    }
    return d || 'M0 0';
  };
  const opacity = (c) => (c.alpha == null ? 1 : c.alpha).toFixed(3);
  const partialsFirst = (curves) => [...curves].sort((a, b) => (a.total ? 1 : 0) - (b.total ? 1 : 0)).filter(visible);
  svg += `<g clip-path="url(#clip)" fill="none" stroke-linejoin="round">`;
  // region under the mouse, and the region whose distances are highlighted in the viewer
  svg += `<rect id="region" x="0" y="${mt}" width="0" height="${H - mt - mb}" fill="#2f6fdb" fill-opacity="0.09" stroke="#2f6fdb" stroke-opacity="0.4" stroke-dasharray="2 3" visibility="hidden"/>`;
  svg += `<rect id="hlRegion" x="0" y="${mt}" width="0" height="${H - mt - mb}" fill="#e69f00" fill-opacity="0.16" stroke="#e69f00" stroke-opacity="0.8" visibility="hidden"/>`;
  if (S.showRef) data.refs.forEach((y) => { svg += `<path d="${path(data.x, y)}" stroke="#d55e00" stroke-width="1.2" stroke-dasharray="5 4"/>`; });
  if (S.showHist) {
    // the total behind its partials, all translucent so that overlapping bars stay visible
    [...data.hist].filter(visible).sort((a, b) => (b.total ? 1 : 0) - (a.total ? 1 : 0)).forEach((c) => {
      svg += `<path d="${bars(data.hx, c.y)}" fill="${c.color}" fill-opacity="${c.total ? 0.22 : 0.55}" stroke="none" opacity="${opacity(c)}"/>`;
    });
  } else {
    partialsFirst(data.curves).forEach((c) => {
      svg += `<path d="${path(data.x, c.y)}" stroke="${c.color}" stroke-width="${c.total ? 2 : 1.5}" opacity="${opacity(c)}"/>`;
    });
  }
  svg += `<line id="xhair" x1="0" x2="0" y1="${mt}" y2="${H - mb}" stroke="#4a5263" stroke-dasharray="3 3" visibility="hidden"/>`;
  svg += `<rect id="zoomBox" x="0" y="${mt}" width="0" height="${H - mt - mb}" fill="rgba(47,111,219,.12)" stroke="#2f6fdb" visibility="hidden"/>`;
  svg += '</g>';
  // legend: click an entry to hide/show its curve (hidden curves are drawn faded in the legend)
  const items = data.curves.filter((c) => c.alpha == null || c.alpha > 0.5)
    .sort((a, b) => (b.total ? 1 : 0) - (a.total ? 1 : 0));
  if (items.length) {
    const lw = 34 + Math.max(...items.map((c) => c.label.length)) * 7.2, lh = 18;
    const lx = conv.legend === 'left' ? ml + 10 : W - mr - lw - 10, ly = mt + 8;
    svg += `<rect x="${lx}" y="${ly}" width="${lw}" height="${items.length * lh + 8}" rx="6" fill="#ffffff" fill-opacity=".9" stroke="#e2e5ea"/>`;
    items.forEach((c, k) => {
      const yy = ly + 4 + k * lh + lh / 2, op = visible(c) ? 1 : 0.35;
      svg += `<g class="lg" data-label="${esc(c.label)}" opacity="${op}"><rect x="${lx + 2}" y="${yy - lh / 2}" width="${lw - 4}" height="${lh}" fill="transparent"/>`;
      svg += `<line x1="${lx + 8}" x2="${lx + 26}" y1="${yy}" y2="${yy}" stroke="${c.color}" stroke-width="${c.total ? 2.5 : 2}"/>`;
      svg += `<text x="${lx + 32}" y="${yy + 4}" fill="${INK}">${c.label === 'total' ? nameSvg(conv.name) : esc(c.label)}</text></g>`;
    });
  }
  svg += '</svg>';
  box.innerHTML = svg;
  drawHlRegion();
}

function redraw() {
  // an animation that is overdue (e.g. paused in a background tab) is finished here
  if (tween && performance.now() - tween.t0 > TWEEN_MS) tween = null;
  if (current && !tween) { shown = current; drawChart(shown); }
}
new ResizeObserver(() => { if (shown) drawChart(shown); }).observe($('chart'));

// ---------- mouse: hover readout, drag zoom, wheel zoom, double-click reset, legend clicks
const chart = $('chart');
const invX = (px) => plot.x0 + ((px - plot.ml) / (plot.W - plot.ml - plot.mr)) * (plot.x1 - plot.x0);
const inPlot = (px, py) => plot && px >= plot.ml && px <= plot.W - plot.mr && py >= plot.mt && py <= plot.H - plot.mb;
let dragFrom = null;
function localXY(ev) { const b = chart.getBoundingClientRect(); return [ev.clientX - b.left, ev.clientY - b.top]; }

// dragging across the plot zooms into the dragged r range (zoom mode, the default) or moves the
// visible range (pan mode); a click without dragging highlights distances in both modes
let panStart = null;
function setDragMode(mode) {
  S.dragMode = mode;
  $('modeZoom').classList.toggle('on', mode === 'zoom');
  $('modePan').classList.toggle('on', mode === 'pan');
  chart.classList.toggle('pan', mode === 'pan');
}
$('modeZoom').onclick = () => setDragMode('zoom');
$('modePan').onclick = () => setDragMode('pan');
chart.addEventListener('mousedown', (ev) => {
  const [px, py] = localXY(ev);
  if (ev.button !== 0 || !inPlot(px, py) || ev.target.closest('.lg')) return;
  ev.preventDefault();
  dragFrom = px;
  const range = xRange(current);
  panStart = { range, perPx: (range[1] - range[0]) / (plot.W - plot.ml - plot.mr) };
  if (S.dragMode === 'pan') chart.classList.add('dragging');
});
window.addEventListener('mousemove', (ev) => {
  if (!plot || !current) return;
  const [px, py] = localXY(ev);
  if (dragFrom != null && S.dragMode === 'pan') {
    // move the range with the mouse, within 0 … r max
    const width = panStart.range[1] - panStart.range[0];
    let a = panStart.range[0] - (px - dragFrom) * panStart.perPx;
    a = Math.max(0, Math.min(current.rMax - width, a));
    if (S.xr) setXRange([a, a + width]);
  } else if (dragFrom != null) {
    const zb = chart.querySelector('#zoomBox');
    const a = Math.max(plot.ml, Math.min(dragFrom, px)), b = Math.min(plot.W - plot.mr, Math.max(dragFrom, px));
    zb.setAttribute('x', a); zb.setAttribute('width', b - a); zb.setAttribute('visibility', 'visible');
  }
  showReadout(px, py);
});
window.addEventListener('mouseup', (ev) => {
  if (dragFrom == null) return;
  const [px, py] = localXY(ev);
  const from = dragFrom;
  dragFrom = null;
  chart.classList.remove('dragging');
  if (Math.abs(px - from) > 5) {
    if (S.dragMode !== 'pan') {
      const a = invX(Math.min(px, from)), b = invX(Math.max(px, from));
      setXRange([Math.max(0, a), Math.min(current.rMax, b)]);
    }
  } else if (inPlot(px, py)) clickHighlight(px, py);
  const zb = chart.querySelector('#zoomBox');
  if (zb) zb.setAttribute('visibility', 'hidden');
});
chart.addEventListener('dblclick', () => setXRange(null));
chart.addEventListener('wheel', (ev) => {
  const [px, py] = localXY(ev);
  if (!current || !inPlot(px, py)) return;
  ev.preventDefault();
  const c = invX(px), f = Math.exp(ev.deltaY * 0.002);
  const [x0, x1] = xRange(current);
  let a = c - (c - x0) * f, b = c + (x1 - c) * f;
  if (b - a < 0.2) return;
  if (a < 0) { b -= a; a = 0; }
  if (b > current.rMax) { a -= b - current.rMax; b = current.rMax; }
  setXRange([Math.max(0, a), b]);
}, { passive: false });
chart.addEventListener('click', (ev) => {
  const g = ev.target.closest('.lg');
  if (g) toggleCurve(g.dataset.label);
});
chart.addEventListener('mouseleave', () => hideReadout());
document.addEventListener('keydown', (ev) => { if (ev.key === 'Escape' && S.hl) clearHighlight(); });

function nearestIndex(xs, r) {
  const i = Math.round((r - xs[0]) / (xs[1] - xs[0]));
  return Math.max(0, Math.min(xs.length - 1, i));
}
function hideReadout() {
  $('readout').hidden = true;
  ['#xhair', '#region'].forEach((id) => { const el = chart.querySelector(id); if (el) el.setAttribute('visibility', 'hidden'); });
}
function setBand(el, a, b) {
  const ra = Math.max(plot.X(a), plot.ml), rb = Math.min(plot.X(b), plot.W - plot.mr);
  el.setAttribute('x', ra); el.setAttribute('width', Math.max(rb - ra, 1));
  el.setAttribute('visibility', rb > ra - 1 ? 'visible' : 'hidden');
}
function inHighlight(r) {
  if (!S.hl) return false;
  const [a, b] = evalRegion(S.hl.r);
  return r >= a && r <= b;
}
// hovering is cheap: values of the curves, the region under the mouse and a hint; the
// distances are only evaluated (and drawn in the viewer) when the user clicks
function showReadout(px, py) {
  if (!inPlot(px, py) || dragFrom != null || tween) { hideReadout(); return; }
  const ro = $('readout'), xh = chart.querySelector('#xhair'), rg = chart.querySelector('#region');
  const d = mainLayer(shown), i = nearestIndex(d.x, invX(px)), x = d.x[i];
  const [a, b] = evalRegion(x);
  const xp = plot.X(x);
  xh.setAttribute('x1', xp); xh.setAttribute('x2', xp); xh.setAttribute('visibility', 'visible');
  setBand(rg, a, b);
  const fmt = (v) => (Math.abs(v) >= 1000 ? v.toFixed(0) : Math.abs(v) >= 1 ? v.toFixed(3) : v.toPrecision(3));
  const rows = d.curves.filter(visible).sort((p, q) => (q.total ? 1 : 0) - (p.total ? 1 : 0))
    .map((c) => `<div><i style="background:${c.color === INK ? '#f3f5f8' : c.color}"></i>${c.total ? nameHtml(CONV[shown.conv].name) : esc(c.label)}: <b>${fmt(c.y[i])}</b></div>`);
  const what = S.showHist ? ' (histogram bin)' : '';
  const hint = inHighlight(x) ? 'left-click: remove the highlighting'
    : `left-click: highlight distances in ${a.toFixed(2)}–${b.toFixed(2)} Å`;
  ro.innerHTML = `<div><i style="background:transparent"></i><b>r = ${x.toFixed(3)} Å</b>${what}</div>${rows.join('')}<div class="hint">${hint}</div>`;
  ro.hidden = false;
  const box = $('pdfWrap').getBoundingClientRect(), cb = chart.getBoundingClientRect();
  const w = ro.offsetWidth, h = ro.offsetHeight;
  let left = xp + 16;
  if (left + w > cb.width - 8) left = xp - w - 16;
  const top = Math.min(Math.max(py - h / 2, 4), cb.height - h - 4) + (cb.top - box.top);
  ro.style.left = `${left}px`; ro.style.top = `${top}px`;
}
// A click highlights the distances of the region at that r, and the highlighting stays while the
// structure is turned; a click into the highlighted region (or Esc) removes it.
function clickHighlight(px, py) {
  const d = mainLayer(shown), r = d.x[nearestIndex(d.x, invX(px))];
  if (inHighlight(r)) clearHighlight();
  else { S.hl = { r }; updateHighlight(); }
  showReadout(px, py);
}
function clearHighlight() {
  S.hl = null;
  drawHlRegion();
  highlightDistances([]);
}
function drawHlRegion() {
  updateControlStates();
  const el = chart.querySelector('#hlRegion');
  if (!el) return;
  if (!S.hl || !S.active) { el.setAttribute('visibility', 'hidden'); return; }
  setBand(el, ...evalRegion(S.hl.r));
}

function setXRange(xr) {
  S.xr = xr && xr[1] - xr[0] > 0.05 ? xr : null;
  syncRange();
  redraw();
  updateControlStates();
}
function syncRange() {
  if (!current) { $('rFrom').value = ''; $('rTo').value = ''; return; }
  const [a, b] = xRange(current);
  $('rFrom').value = +a.toFixed(2); $('rTo').value = +b.toFixed(2);
  $('rFrom').max = $('rTo').max = current.rMax;
}
$('rFrom').onchange = $('rTo').onchange = () => {
  if (!current) return;
  const a = Math.max(0, +$('rFrom').value), b = Math.min(current.rMax, +$('rTo').value);
  setXRange(b > a ? [a, b] : null);
};
$('resetRange').onclick = $('resetZoom').onclick = () => setXRange(null);

// ------------------------------------------------------------------ PDF options
document.querySelectorAll('[data-conv]').forEach((b) => {
  b.onclick = () => {
    if (S.conv === b.dataset.conv && S.rad === b.dataset.rad) return;
    const radChanged = S.rad !== b.dataset.rad;
    S.conv = b.dataset.conv; S.rad = b.dataset.rad;
    updatePdf(true);
    if (radChanged && S.hl) updateHighlight();    // the shares depend on the scattering lengths
  };
});
function setSigma(v, fromRange) {
  if (!(v > 0)) return;
  S.sigma = v;
  if (fromRange) $('sigma').value = v; else $('sigmaRange').value = v;
  updatePdf(false);
  if (S.hl) updateHighlight();
}
$('sigmaRange').oninput = () => setSigma(+$('sigmaRange').value, true);
$('sigma').onchange = () => setSigma(+$('sigma').value, false);
// broadened function and histogram are exclusive (two radio buttons)
function setLayer(hist) {
  S.showHist = hist;
  $('histogram').checked = hist;
  $('showBroad').checked = !hist;
  redraw();
  if (S.hl) updateHighlight();
}
$('showBroad').onchange = () => setLayer(false);
$('histogram').onchange = () => setLayer(true);

function toggleCurve(label) {
  if (S.hidden.has(label)) S.hidden.delete(label); else S.hidden.add(label);
  renderCurveList();
  redraw();
  if (S.hl) updateHighlight();
}
function partialLabels() { return S.active ? Object.keys(S.active.pdf.partials) : []; }
$('curvesAll').onclick = () => { S.hidden.clear(); renderCurveList(); redraw(); if (S.hl) updateHighlight(); };
$('curvesNone').onclick = () => { partialLabels().forEach((l) => S.hidden.add(l)); renderCurveList(); redraw(); if (S.hl) updateHighlight(); };
// highlighting options: half-width of the region and minimum share of the R(r) area
$('hlWidth').onchange = () => {
  const v = +$('hlWidth').value;
  if (v > 0) S.hlWidth = v; else $('hlWidth').value = S.hlWidth;
  redraw();
  if (S.hl) updateHighlight();
};
$('hlMin').onchange = () => {
  const v = +$('hlMin').value;
  if (v >= 0 && v <= 100) S.hlMin = v / 100; else $('hlMin').value = S.hlMin * 100;
  if (S.hl) updateHighlight();
};
$('hlClear').onclick = () => clearHighlight();
// appearance of the highlighting: rod thickness and visibility of the other atoms
function showHlLook() {
  $('hlRodVal').textContent = `${S.hlRod.toFixed(2)} Å`;
  $('hlDimVal').textContent = S.hlDim > 0 ? `${Math.round(S.hlDim * 100)} %` : 'hidden';
}
$('hlRod').oninput = () => {
  S.hlRod = +$('hlRod').value;
  showHlLook();
  if (S.hl) { hlKey = ''; updateHighlight(); }
};
$('hlDim').oninput = () => { S.hlDim = +$('hlDim').value; showHlLook(); applyDim(); };
showHlLook();

function renderCurveList() {
  const list = $('curveList');
  if (!S.active) { list.className = 'muted'; list.textContent = '–'; updateControlStates(); return; }
  list.className = '';
  const conv = CONV[S.conv];
  const row = (checked, swatch, name, onchange, { off = false, title = '' } = {}) => {
    const lab = document.createElement('label');
    lab.className = `${checked ? '' : 'unticked'}${off ? ' off' : ''}`;
    lab.title = title;
    lab.innerHTML = `<input type="checkbox" ${checked ? 'checked' : ''} ${off ? 'disabled' : ''}>${swatch}<span>${name}</span>`;
    lab.querySelector('input').onchange = onchange;
    return lab;
  };
  const sw = (color, cls = '') => `<i class="sw ${cls}" style="background:${color}"></i>`;
  const partialTitle = conv.totalOnly ? `${conv.name} is not weighted by scattering lengths and is shown as total only` : '';
  list.replaceChildren(
    row(!S.hidden.has('total'), sw(INK, 'total'), `${nameHtml(conv.name)} (total)`, () => toggleCurve('total')),
    ...partialLabels().map((l, i) => row(!S.hidden.has(l), sw(PALETTE[i % PALETTE.length]), esc(l), () => toggleCurve(l),
      { off: !!conv.totalOnly, title: partialTitle })),
    row(S.showRef, '<i class="sw dash"></i>', 'asymptotes', () => { S.showRef = !S.showRef; renderCurveList(); redraw(); },
      { title: 'The low-r and high-r limits of the selected function (Peterson et al. 2021, Table 4), drawn as dashed lines, e.g. G(r) → −4πρ₀r below the first peak and → 0 at high r. Useful to check the normalisation.' }));
  updateControlStates();
}

// ------------------------------------------------------------------ export
function baseName() {
  const e = S.active;
  const file = e.file.replace(/\.cif$/i, '');
  return (e.name === file ? file : `${file}_${e.name}`).replace(/[^\w.-]+/g, '_');
}
function fileKey() { return `${CONV[current.conv].file}${current.rad === 'neutron' ? '_neutron' : ''}`; }
function download(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
$('exportCsv').onclick = () => {
  if (!current) return;
  // the curves on screen: the broadened function or the histogram
  const { x: xs, curves } = mainLayer(current);
  const header = ['radii', ...curves.map((c) => c.label)].join(',');
  const lines = Array.from(xs, (x, i) => [+x.toFixed(6), ...curves.map((c) => +c.y[i].toPrecision(7))].join(','));
  download(new Blob([`${header}\n${lines.join('\n')}\n`], { type: 'text/csv' }), `${baseName()}_${fileKey()}.csv`);
};
function chartSvgText() {
  const svg = $('chart').querySelector('svg').cloneNode(true);
  svg.querySelectorAll('#xhair, #zoomBox, #region, #hlRegion').forEach((el) => el.remove());
  return new XMLSerializer().serializeToString(svg);
}
$('exportSvg').onclick = () => {
  if (!current) return;
  download(new Blob([chartSvgText()], { type: 'image/svg+xml' }), `${baseName()}_${fileKey()}.svg`);
};
$('exportPng').onclick = () => {
  if (!current) return;
  const img = new Image(), scale = 3;
  img.onload = () => {
    const cv = document.createElement('canvas');
    cv.width = plot.W * scale; cv.height = plot.H * scale;
    const g = cv.getContext('2d');
    g.scale(scale, scale); g.drawImage(img, 0, 0, plot.W, plot.H);
    cv.toBlob((b) => download(b, `${baseName()}_${fileKey()}.png`));
  };
  img.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(chartSvgText())}`;
};

// ------------------------------------------------------------------ splitter between viewer and PDF
$('splitter').addEventListener('pointerdown', (ev) => {
  const center = $('center'), sp = $('splitter');
  sp.setPointerCapture(ev.pointerId); sp.classList.add('drag');
  const move = (e) => {
    const b = center.getBoundingClientRect();
    const h = Math.min(b.height - 120, Math.max(160, b.bottom - e.clientY));
    center.style.setProperty('--pdf-h', `${h}px`);
  };
  const up = () => { sp.classList.remove('drag'); sp.removeEventListener('pointermove', move); sp.removeEventListener('pointerup', up); };
  sp.addEventListener('pointermove', move); sp.addEventListener('pointerup', up);
});

// ------------------------------------------------------------------ three.js viewer
const viewer = $('viewer');
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
viewer.appendChild(renderer.domElement);
// the scene is drawn only when something has changed, so an idle viewer costs no GPU time
let needsRender = true;
const requestRender = () => { needsRender = true; };
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(30, 1, 0.1, 50000);
camera.position.set(30, 25, 40);
const controls = new OrbitControls(camera, renderer.domElement);
controls.addEventListener('change', requestRender);
scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const keyLight = new THREE.DirectionalLight(0xffffff, 1.6);
keyLight.position.set(5, 8, 10);
camera.add(keyLight);
scene.add(camera);
const root = new THREE.Group();
scene.add(root);

new ResizeObserver(() => {
  const w = viewer.clientWidth, h = viewer.clientHeight;
  renderer.setSize(w, h);
  camera.aspect = w / Math.max(h, 1);
  camera.updateProjectionMatrix();
  requestRender();
}).observe(viewer);

// crystallographic axis compass (bottom right), drawn in a colour that contrasts with the background
const axScene = new THREE.Scene();
const axCam = new THREE.OrthographicCamera(-1.45, 1.45, 1.45, -1.45, -10, 10);
const axArrows = [0, 1, 2].map(() => {
  const arr = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 1, '#ffffff', 0.3, 0.17);
  axScene.add(arr);
  return arr;
});
function labelTexture(text, col) {
  const cv = document.createElement('canvas'); cv.width = cv.height = 128;
  const g = cv.getContext('2d');
  g.font = 'italic bold 92px Georgia, "Times New Roman", serif'; g.textAlign = 'center'; g.textBaseline = 'middle';
  g.fillStyle = col; g.fillText(text, 64, 70);
  return new THREE.CanvasTexture(cv);
}
const axLabels = ['a', 'b', 'c'].map(() => {
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ depthTest: false }));
  sp.scale.set(0.62, 0.62, 1);
  axScene.add(sp);
  return sp;
});
function setBackground() {
  const white = $('whiteBg').checked, ink = white ? '#000000' : '#ffffff';
  scene.background = new THREE.Color(white ? '#ffffff' : '#11151d');
  $('viewerWrap').style.background = white ? '#ffffff' : 'var(--viewer)';
  axArrows.forEach((a) => a.setColor(ink));
  axLabels.forEach((sp, i) => { sp.material.map?.dispose(); sp.material.map = labelTexture('abc'[i], ink); sp.material.needsUpdate = true; });
  if (cellLines) cellLines.material.color.set(white ? '#444' : '#9aa5b8');
  requestRender();
}
function renderCompass() {
  const size = 130, w = viewer.clientWidth;
  axCam.position.copy(camera.position).sub(controls.target).normalize().multiplyScalar(5);
  axCam.up.copy(camera.up);
  axCam.lookAt(0, 0, 0);
  renderer.autoClear = false;
  renderer.clearDepth();
  renderer.setScissorTest(true);
  renderer.setScissor(w - size - 6, 6, size, size);
  renderer.setViewport(w - size - 6, 6, size, size);
  renderer.render(axScene, axCam);
  renderer.setScissorTest(false);
  renderer.setViewport(0, 0, w, viewer.clientHeight);
  renderer.autoClear = true;
}

const sphereGeo = new THREE.SphereGeometry(1, 20, 14);
const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 10, 1);
const solidMat = new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0.05 });
// bonds to partly occupied sites are translucent
const partialMat = new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0.05, transparent: true, opacity: 0.45, depthWrite: false });
// Atoms: the occupancy of a site is shown like a cake. The sphere is coloured over the angle
// around the vertical axis up to occupancy x 360°; the sphere's area is spread evenly over that
// angle, so the coloured fraction of the surface equals the occupancy. On mixed sites the second
// element follows as a slice in its own colour. The empty slice is cut out, and two half-disc
// walls (wallMat, drawn only for partly occupied sites) close the cut.
const TAU = 6.283185307179586;
const atomMat = new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0.05 });
atomMat.onBeforeCompile = (shader) => {
  shader.vertexShader = shader.vertexShader
    .replace('#include <common>', `#include <common>
attribute float occA;
attribute float occT;
attribute vec3 colB;
varying vec3 vLocal;
varying float vOccA;
varying float vOccT;
varying vec3 vColB;`)
    .replace('#include <begin_vertex>', `#include <begin_vertex>
vLocal = position;
vOccA = occA; vOccT = occT; vColB = colB;`);
  shader.fragmentShader = shader.fragmentShader
    .replace('#include <common>', `#include <common>
varying vec3 vLocal;
varying float vOccA;
varying float vOccT;
varying vec3 vColB;`)
    .replace('#include <color_fragment>', `#include <color_fragment>
float turn = fract(atan(vLocal.z, vLocal.x) / ${TAU} + 1.0);
if (turn >= vOccT) discard;
if (turn >= vOccA) diffuseColor.rgb = vColB;`);
};
// the two walls of the cut: a half disc (x >= 0 in the xy plane) turned about the vertical axis
// to the angle 0 (wall 0) or occupancy x 360° (wall 1); a slightly darker 'inside' colour
const wallMat = new THREE.MeshStandardMaterial({ roughness: 0.6, metalness: 0.0, side: THREE.DoubleSide });
wallMat.onBeforeCompile = (shader) => {
  shader.vertexShader = shader.vertexShader
    .replace('#include <common>', `#include <common>
attribute float wall;
attribute float occT;
attribute vec3 colB;
varying float vWall;
varying vec3 vColB;
float wallAngle;`)
    .replace('#include <beginnormal_vertex>', `#include <beginnormal_vertex>
wallAngle = wall > 0.5 ? occT * ${TAU} : 0.0;
objectNormal = vec3(-sin(wallAngle), 0.0, cos(wallAngle));`)
    .replace('#include <begin_vertex>', `#include <begin_vertex>
transformed = vec3(position.x * cos(wallAngle), position.y, position.x * sin(wallAngle));
vWall = wall; vColB = colB;`);
  shader.fragmentShader = shader.fragmentShader
    .replace('#include <common>', `#include <common>
varying float vWall;
varying vec3 vColB;`)
    .replace('#include <color_fragment>', `#include <color_fragment>
if (vWall > 0.5) diffuseColor.rgb = vColB;
diffuseColor.rgb *= 0.8;`);
};
// geometry of the walls: two copies of a half disc, marked by the attribute 'wall'
const wallGeo = (() => {
  const half = new THREE.CircleGeometry(1, 24, -Math.PI / 2, Math.PI);
  const pos = half.attributes.position.array, idx = half.index.array, n = pos.length / 3;
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array([...pos, ...pos]), 3));
  geo.setAttribute('normal', new THREE.BufferAttribute(new Float32Array(6 * n).fill(0), 3));
  geo.setAttribute('wall', new THREE.BufferAttribute(new Float32Array(2 * n).fill(1, n), 1));
  geo.setIndex([...idx, ...Array.from(idx, (i) => i + n)]);
  half.dispose();
  return geo;
})();
let meshes = [], cellLines = null;
const MAX_ATOMS = 40000;

const vadd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const vscale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const imgVec = (n, L) => [0, 1, 2].map((k) => n[0] * L[0][k] + n[1] * L[1][k] + n[2] * L[2][k]);
const _m = new THREE.Matrix4(), _q = new THREE.Quaternion(), _p = new THREE.Vector3(), _s = new THREE.Vector3();
const _Y = new THREE.Vector3(0, 1, 0), _d = new THREE.Vector3(), _c = new THREE.Color();

function atomRadius(a) { return a.el === 'H' ? 0.17 : 0.2 + 0.22 * a.r; }
function atomShown(a) { return a.el !== 'H' && a.el !== 'D' || $('showH').checked; }

function instanced(geo, mat, items, place) {
  if (!items.length) return null;
  const mesh = new THREE.InstancedMesh(geo, mat, items.length);
  items.forEach((it, k) => { place(it); mesh.setMatrixAt(k, _m.compose(_p, _q, _s)); mesh.setColorAt(k, _c.set(it.color)); });
  mesh.userData.items = items;
  mesh.computeBoundingSphere();
  root.add(mesh);
  meshes.push(mesh);
  return mesh;
}

function buildScene() {
  highlightDistances([]);
  meshes.forEach((m) => { root.remove(m); if (m.userData.atoms || m.userData.ownGeometry) m.geometry.dispose(); m.dispose(); });
  meshes = [];
  if (cellLines) { root.remove(cellLines); cellLines.geometry.dispose(); cellLines = null; }
  requestRender();
  renderLegend();
  const e = S.active;
  if (!e) return;
  const st = e.st, L = st.L;
  const imgs = [];
  for (let i = 0; i < S.sc[0]; i++) for (let j = 0; j < S.sc[1]; j++) for (let k = 0; k < S.sc[2]; k++) imgs.push([i, j, k]);
  const atoms = [], bonds = [];
  const inside = (m) => m.every((v, k) => v >= 0 && v < S.sc[k]);
  const boundary = $('boundaryAtoms').checked, showBonds = $('showBonds').checked;
  const ghosts = new Set();
  const addAtom = (s, m) => {
    const a = st.atoms[s];
    atoms.push({ s, x: vadd(a.xyz, imgVec(m, L)), color: a.color });
  };
  const halfBond = (s, ms, t, mt) => {
    const as = st.atoms[s], at = st.atoms[t];
    const xs = vadd(as.xyz, imgVec(ms, L)), xt = vadd(at.xyz, imgVec(mt, L));
    bonds.push({ a: xs, b: vscale(vadd(xs, xt), 0.5), color: as.color, partial: as.occ < 0.99 || at.occ < 0.99 });
  };
  for (const img of imgs) {
    st.atoms.forEach((a, s) => { if (atomShown(a)) addAtom(s, img); });
    // every bond is listed in both directions, so each atom inside the cells draws the half next to
    // it; a partner outside the cells is drawn as a boundary atom (once) with its half of the bond,
    // or, without boundary atoms, the bond is left out, so that no bond ends in empty space
    st.bonds.forEach(([i, j, n]) => {
      if (!atomShown(st.atoms[i]) || !atomShown(st.atoms[j])) return;
      const m = vadd(img, n);
      if (inside(m)) { if (showBonds) halfBond(i, img, j, m); return; }
      if (!boundary) return;
      const key = `${j}:${m}`;
      if (!ghosts.has(key)) { ghosts.add(key); addAtom(j, m); }
      if (showBonds) { halfBond(i, img, j, m); halfBond(j, m, i, img); }
    });
  }
  if (atoms.length > MAX_ATOMS * 1.5) {
    // boundary atoms can push a large supercell over the limit: fall back to one cell
    showSupercell(`<span class="err">${atoms.length.toLocaleString()} atoms with boundary atoms: showing one cell</span>`);
    S.sc.splice(0, 3, 1, 1, 1);
    buildScene();
    return;
  }
  const placeAtom = (it) => { const r = atomRadius(st.atoms[it.s]); _p.set(...it.x); _q.identity(); _s.set(r, r, r); };
  const placeBond = (it) => {
    _d.set(it.b[0] - it.a[0], it.b[1] - it.a[1], it.b[2] - it.a[2]);
    const len = _d.length() || 1e-6;
    _q.setFromUnitVectors(_Y, _d.multiplyScalar(1 / len));
    _p.set((it.a[0] + it.b[0]) / 2, (it.a[1] + it.b[1]) / 2, (it.a[2] + it.b[2]) / 2);
    _s.set(0.11, len, 0.11);
  };
  // atoms: one mesh with its own copy of the sphere geometry, which carries the occupancies
  const atomMesh = instanced(sphereGeo.clone(), atomMat, atoms, placeAtom);
  if (atomMesh) {
    const n = atoms.length, occA = new Float32Array(n), occT = new Float32Array(n), colB = new Float32Array(3 * n);
    atoms.forEach((it, k) => {
      const parts = st.atoms[it.s].parts || [[st.atoms[it.s].el, st.atoms[it.s].occ, it.color]];
      occA[k] = parts[0][1];
      occT[k] = Math.min(1, parts.reduce((sum, p) => sum + p[1], 0));
      _c.set(parts[1] ? parts[1][2] : parts[0][2]);
      colB.set([_c.r, _c.g, _c.b], 3 * k);
      if (occT[k] > 0.995) occT[k] = 1.01;     // fully occupied: nothing cut out
      if (!parts[1] && occA[k] > 0.995) occA[k] = 1.01;
    });
    atomMesh.geometry.setAttribute('occA', new THREE.InstancedBufferAttribute(occA, 1));
    atomMesh.geometry.setAttribute('occT', new THREE.InstancedBufferAttribute(occT, 1));
    atomMesh.geometry.setAttribute('colB', new THREE.InstancedBufferAttribute(colB, 3));
    atomMesh.userData.atoms = true;
    // walls of the cut for the partly occupied sites (own geometry copy for the attributes)
    const cut = atoms.map((it, k) => k).filter((k) => occT[k] < 1);
    const wallMesh = instanced(wallGeo.clone(), wallMat, cut.map((k) => atoms[k]), placeAtom);
    if (wallMesh) {
      wallMesh.geometry.setAttribute('occT', new THREE.InstancedBufferAttribute(Float32Array.from(cut, (k) => occT[k]), 1));
      wallMesh.geometry.setAttribute('colB', new THREE.InstancedBufferAttribute(Float32Array.from(cut.flatMap((k) => [...colB.subarray(3 * k, 3 * k + 3)])), 3));
      wallMesh.userData.ownGeometry = true;
    }
  }
  instanced(cylGeo, solidMat, bonds.filter((b) => !b.partial), placeBond);
  instanced(cylGeo, partialMat, bonds.filter((b) => b.partial), placeBond);
  // unit cell box of the whole supercell
  if ($('showCell').checked) {
    const corners = [];
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) for (let k = 0; k < 2; k++) corners.push(imgVec([i * S.sc[0], j * S.sc[1], k * S.sc[2]], L));
    const edges = [[0, 1], [2, 3], [4, 5], [6, 7], [0, 2], [1, 3], [4, 6], [5, 7], [0, 4], [1, 5], [2, 6], [3, 7]];
    const pos = new Float32Array(72);
    edges.forEach(([a, b], k) => { pos.set(corners[a], k * 6); pos.set(corners[b], k * 6 + 3); });
    cellLines = new THREE.LineSegments(new THREE.BufferGeometry(),
      new THREE.LineBasicMaterial({ color: $('whiteBg').checked ? '#444' : '#9aa5b8' }));
    cellLines.geometry.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    root.add(cellLines);
  }
  L.forEach((v, i) => {
    const d = new THREE.Vector3(...v).normalize();
    axArrows[i].setDirection(d);
    axLabels[i].position.copy(d).multiplyScalar(1.28);
  });
  if (S.hl) updateHighlight();             // the highlighting stays when the scene changes
}

function renderLegend() {
  const e = S.active;
  if (!e) { $('legend').innerHTML = ''; return; }
  const seen = new Map();
  // every element, also those that are only the minor part of mixed sites
  e.st.atoms.forEach((a) => (a.parts || [[a.el, a.occ, a.color]]).forEach(([el, , color]) => { if (!seen.has(el)) seen.set(el, color); }));
  let html = [...seen].map(([el, c]) => `<i style="background:${c}"></i>${esc(el)}`).join('');
  if (e.st.atoms.some((a) => a.occ < 0.99 || (a.parts && a.parts.length > 1)))
    html += '<i class="occ"></i>partly occupied (coloured slice = occupancy)';
  $('legend').innerHTML = html;
}

function frameCamera(dir) {
  const e = S.active;
  if (!e) return;
  const L = e.st.L;
  const center = vscale(imgVec(S.sc, L), 0.5);
  const diag = Math.hypot(...imgVec(S.sc, L));
  const d = new THREE.Vector3(...(dir || [0.6, 0.45, 1])).normalize();
  controls.target.set(...center);
  camera.up.set(0, 1, 0);
  if (Math.abs(d.y) > 0.95) camera.up.set(0, 0, 1);
  const dist = (1.1 * diag) / (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov / 2)));
  camera.position.copy(controls.target).addScaledVector(d, dist);
  camera.near = Math.max(0.05, dist / 200); camera.far = dist + 50 * diag;
  camera.updateProjectionMatrix();
  controls.update();
  requestRender();
}
document.querySelectorAll('[data-view]').forEach((b) => {
  b.onclick = () => {
    const e = S.active;
    if (!e) return;
    const L = e.st.L, k = 'abc'.indexOf(b.dataset.view);
    // look along the lattice vector (down the axis towards the origin)
    frameCamera(L[k]);
  };
});
// ---------- supercell: a x b x c repeats, limited by the number of drawn atoms
const MAX_SC = 12;
function visibleAtomCount() { return S.active ? S.active.st.atoms.filter(atomShown).length : 0; }
function showSupercell(note) {
  ['scA', 'scB', 'scC'].forEach((id, k) => { $(id).value = S.sc[k]; });
  const n = visibleAtomCount() * S.sc[0] * S.sc[1] * S.sc[2];
  $('scInfo').innerHTML = note || (S.active ? `${n.toLocaleString()} atoms` : '');
}
function setSupercell(sc) {
  sc = sc.map((v) => Math.min(MAX_SC, Math.max(1, Math.round(v) || 1)));
  const n = visibleAtomCount() * sc[0] * sc[1] * sc[2];
  if (n > MAX_ATOMS) {
    showSupercell(`<span class="err">${sc.join('×')} = ${n.toLocaleString()} atoms (max ${MAX_ATOMS.toLocaleString()})</span>`);
    return false;
  }
  const changed = sc.some((v, k) => v !== S.sc[k]);
  S.sc.splice(0, 3, ...sc);               // the array belongs to the structure on screen
  showSupercell();
  if (changed) { buildScene(); frameCamera(); }
  return true;
}
// thickness of one cell along each lattice direction: the spacing of the lattice planes spanned
// by the other two vectors, 1/|b_k| (b: reciprocal vectors); correct for triclinic cells too
function cellThickness(L) {
  const Linv = inverse3(L);
  return [0, 1, 2].map((c) => 1 / Math.hypot(Linv[0][c], Linv[1][c], Linv[2][c]));
}
// repeats giving about the same thickness t in every direction; too many atoms: the largest
// repeat is reduced until the supercell fits
function fitSupercell(t) {
  if (!S.active) return;
  let sc = cellThickness(S.active.st.L).map((h) => Math.max(1, Math.min(MAX_SC, Math.round(t / h))));
  while (!setSupercell(sc) && sc.some((v) => v > 1)) {
    const k = sc.indexOf(Math.max(...sc));
    sc = sc.map((v, i) => (i === k ? v - 1 : v));
  }
}
['scA', 'scB', 'scC'].forEach((id) => {
  $(id).onchange = () => setSupercell(['scA', 'scB', 'scC'].map((x) => +$(x).value));
});
$('scCube').onclick = () => { if (S.active) fitSupercell(Math.max(...cellThickness(S.active.st.L))); };
$('scRmax').onclick = () => { if (S.active) fitSupercell(S.active.pdf.r_max); };
['showBonds', 'showCell', 'boundaryAtoms'].forEach((id) => { $(id).onchange = buildScene; });
$('showH').onchange = () => { buildScene(); showSupercell(); };
$('whiteBg').onchange = setBackground;
$('savePng').onclick = () => {
  if (!S.active) return;
  renderer.render(scene, camera);
  renderer.domElement.toBlob((b) => download(b, `${baseName()}_structure.png`));
};

// ---------- distances of a region of the PDF, highlighted in the viewer on click
// The server groups all atom pairs into distance types (symmetry-equivalent pairs, one
// representative each). A click ranks the types in the region by their share of the area under
// R(r) and shows those above the threshold as coloured pairs, with the rest of the structure dimmed.
// Colour-blind safe hues (Okabe & Ito), ordered from bright to dark; the colour of a distance is
// also muted towards grey the smaller its share is compared with the strongest one, so the
// dominant distances stand out
const HL_COLORS = ['#F0E442', '#56B4E9', '#E69F00', '#009E73', '#CC79A7', '#0072B2', '#D55E00', '#BBBBBB'];
const MAX_HIGHLIGHTS = HL_COLORS.length;
const HL_MUTE_TO = '#5b6170', HL_MAX_MUTE = 0.6;
function highlightColor(rank, relative) {
  const mix = HL_MAX_MUTE * (1 - Math.min(1, Math.max(0, relative)));
  return `#${new THREE.Color(HL_COLORS[rank]).lerp(new THREE.Color(HL_MUTE_TO), mix).getHexString()}`;
}
// the distance list is fetched on the first click into a structure's plot
async function ensureDistances(e) {
  if (e.dist) return e.dist;
  if (!e.distRequest) {
    e.distRequest = api('/api/distances', { text: e.text, block: e.block, r_max: e.pdf.r_max })
      .then(({ distances }) => { e.dist = distances; return distances; })
      .finally(() => { e.distRequest = null; });
  }
  return e.distRequest;
}
// region around r: ± the half-width set in box 3; for the histogram, the bins whose centres lie
// in it (at least the bin at r)
function evalRegion(r) {
  const h = S.hlWidth, pdf = S.active.pdf;
  if (!S.showHist) return [r - h, r + h];
  const bin = (x) => Math.round((x - pdf.r0) / pdf.dr);
  const lo = Math.min(bin(r - h + 1e-9), bin(r)), hi = Math.max(bin(r + h - 1e-9), bin(r));
  return [pdf.r0 + (lo - 0.5) * pdf.dr, pdf.r0 + (hi + 0.5) * pdf.dr];
}
// error function (Abramowitz & Stegun 7.1.26, error < 1.5e-7)
function erf(x) {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return x >= 0 ? y : -y;
}
// fraction of the (broadened) contribution of a distance d that falls into [a, b]
function fractionIn(d, a, b) {
  if (S.showHist) return d >= a && d < b ? 1 : 0;
  const z = Math.SQRT2 * S.sigma;
  return 0.5 * (erf((b - d) / z) - erf((a - d) / z));
}
// The distance types in [a, b] whose share of the area under R(r) there reaches the threshold.
// R(r) is used because it is never negative, so the shares of all distances add up to 100 %
// whichever convention is plotted; the area of a type in R(r) is norm * w (see distance_types).
function dominantDistances(a, b) {
  const e = S.active, D = e && e.dist;
  if (!D) return { types: [], others: 0 };
  // total area under R(r) in the region, from the histograms
  // neutron partials can be negative (negative scattering lengths): the shares then refer to
  // the area under |R(r)| summed over the pairs, and the sign of each contribution is shown
  const pdf = e.pdf, neutron = S.rad === 'neutron';
  const W = neutron ? D.wn : D.w, norm = neutron ? D.norm_n : D.norm;
  if (!W || !norm) return { types: [], others: 0 };
  e.Rtot = e.Rtot || {};
  if (!e.Rtot[S.rad]) {
    const { scale } = radiationWeights(pdf, S.rad);
    const sum = new Float64Array(pdf.n);
    Object.entries(pdf.partials).forEach(([l, y]) => { const f = Math.abs(scale(l)); y.forEach((v, k) => { sum[k] += v * f; }); });
    e.Rtot[S.rad] = sum;
  }
  const Rabs = e.Rtot[S.rad];
  const reach = S.showHist ? 0 : 6 * S.sigma;
  let total = 0;
  const k0 = Math.max(0, Math.floor((a - reach - pdf.r0) / pdf.dr)), k1 = Math.min(pdf.n - 1, Math.ceil((b + reach - pdf.r0) / pdf.dr));
  for (let k = k0; k <= k1; k++) total += Rabs[k] * pdf.dr * fractionIn(pdf.r0 + k * pdf.dr, a, b);
  if (!(total > 0)) return { types: [], others: 0 };
  // types are sorted by d: start at the first one that can reach the region
  let k = 0, hi = D.d.length;
  while (k < hi) { const m = (k + hi) >> 1; if (D.d[m] < a - reach) k = m + 1; else hi = m; }
  const cand = [];
  for (; k < D.d.length && D.d[k] <= b + reach; k++) {
    const pair = [e.st.atoms[D.i[k]].el, e.st.atoms[D.j[k]].el].sort().join('-');
    if (S.hidden.has(pair)) continue;     // partials hidden in the plot are not highlighted
    const share = norm * Math.abs(W[k]) * fractionIn(D.d[k], a, b) / total;
    if (share > 0 && share >= S.hlMin) cand.push({ k, pair: pair.replace('-', '–'), d: D.d[k], share, negative: W[k] < 0 });
  }
  cand.sort((p, q) => q.share - p.share);
  const types = cand.slice(0, MAX_HIGHLIGHTS)
    .map((c, rank) => ({ ...c, n: D.n[c.k], color: highlightColor(rank, c.share / cand[0].share) }));
  return { types, others: Math.max(0, 1 - types.reduce((sum, t) => sum + t.share, 0)), more: cand.length - types.length };
}
// (re)evaluates the highlighted region, e.g. after a click or a change of σ or the options
async function updateHighlight() {
  drawHlRegion();
  const e = S.active;
  if (!e || !S.hl) { highlightDistances([]); return; }
  if (!e.dist) {
    $('hlLegend').innerHTML = '<div><span>finding the distances…</span></div>';
    $('hlLegend').style.top = `${$('toolbar').offsetTop + $('toolbar').offsetHeight + 6}px`;
    $('hlLegend').hidden = false;
    try { await ensureDistances(e); } catch (err) { fail(err); }
    if (e !== S.active || !S.hl) return;
  }
  const [a, b] = evalRegion(S.hl.r);
  const result = dominantDistances(a, b);
  highlightDistances(result.types, { a, b, ...result });
}

let hlGroup = null, hlKey = '';
const MAX_HIGHLIGHT_PAIRS = 10000;
// rod radius of the highlighted pairs (S.hlRod, Å) and opacity of the rest of the structure while
// distances are highlighted (S.hlDim, 0: hidden) are set in box 3
let dimmed = false;
function applyDim() {
  const op = dimmed ? S.hlDim : 1;
  [solidMat, atomMat, wallMat].forEach((m) => {
    m.transparent = op < 1; m.opacity = op; m.depthWrite = op >= 1; m.visible = op > 0; m.needsUpdate = true;
  });
  partialMat.opacity = dimmed ? S.hlDim * 0.6 : 0.45;
  partialMat.visible = op > 0;
  partialMat.needsUpdate = true;
  requestRender();
}
function inverse3(m) {
  const [[a, b, c], [d, e, f], [g, h, i]] = m;
  const A = e * i - f * h, B = f * g - d * i, C = d * h - e * g, det = a * A + b * B + c * C;
  return [[A, c * h - b * i, b * f - c * e], [B, a * i - c * g, c * d - a * f], [C, b * g - a * h, a * e - b * d]]
    .map((row) => row.map((v) => v / det));
}
// all pairs of a distance type in the cell: same two site labels at the same distance as its
// representative (symmetry-equivalent pairs), as [i, j, cell offset of j]; cached per structure
function typePairs(e, k) {
  e.pairCache = e.pairCache || new Map();
  if (e.pairCache.has(k)) return e.pairCache.get(k);
  const st = e.st, L = st.L, D = e.dist, d = D.d[k];
  const la = st.atoms[D.i[k]].label, lb = st.atoms[D.j[k]].label;
  const A = [], B = [];
  st.atoms.forEach((a, s) => { if (a.label === la) A.push(s); if (a.label === lb) B.push(s); });
  // cell offsets needed: distance / spacing of the lattice planes, rounded up, + 1
  const Linv = inverse3(L);
  const reach = [0, 1, 2].map((c) => Math.ceil(d * Math.hypot(Linv[0][c], Linv[1][c], Linv[2][c])) + 1);
  const pairs = [];
  for (const i of A) {
    const xi = st.atoms[i].xyz;
    for (const j of B) {
      const xj = st.atoms[j].xyz;
      for (let n0 = -reach[0]; n0 <= reach[0]; n0++) for (let n1 = -reach[1]; n1 <= reach[1]; n1++) for (let n2 = -reach[2]; n2 <= reach[2]; n2++) {
        const t = imgVec([n0, n1, n2], L);
        const dist = Math.hypot(xj[0] + t[0] - xi[0], xj[1] + t[1] - xi[1], xj[2] + t[2] - xi[2]);
        if (Math.abs(dist - d) > 6e-4) continue;
        // like pairs are found from both ends: keep each once
        if (la === lb && (j < i || (j === i && n0 * 1e6 + n1 * 1e3 + n2 < 0))) continue;
        pairs.push([i, j, [n0, n1, n2]]);
      }
    }
  }
  e.pairCache.set(k, pairs);
  return pairs;
}
// draws the pairs of the given distance types in all shown cells; info: region and remainder
function highlightDistances(types, info) {
  const key = S.active && info ? `${S.active.id}:${S.rad}:${types.map((t) => `${t.k}/${t.share.toFixed(4)}`).join(',')}:${info.a}:${info.b}` : '';
  if (key === hlKey) return;
  hlKey = key;
  if (hlGroup) {
    root.remove(hlGroup);
    hlGroup.children.forEach((o) => { o.material.dispose(); o.dispose?.(); });
    hlGroup = null;
  }
  const on = types.length > 0 && !!S.active;
  // dim (or hide) everything else
  dimmed = on;
  applyDim();
  $('hlLegend').hidden = !info;
  requestRender();
  if (!info) return;
  // below the toolbar, which wraps to several lines in narrow windows
  $('hlLegend').style.top = `${$('toolbar').offsetTop + $('toolbar').offsetHeight + 6}px`;
  const pct = (v) => `${(100 * v).toFixed(1)} %`;
  const head = `<div class="lock">${info.a.toFixed(2)}–${info.b.toFixed(2)} Å: share of the ${S.rad === 'neutron' ? 'neutron |R(r)|' : 'X-ray R(r)'} area</div>`;
  const foot = (info.others > 0.0005 ? `<div><i style="background:transparent"></i>other distances ${pct(info.others)}${info.more ? ` (${info.more} more above ${pct(S.hlMin)})` : ''}</div>` : '')
    + '<div><span>click the highlighted region or press Esc to remove</span></div>';
  if (!on) {
    $('hlLegend').innerHTML = `${head}<div><span>no distance reaches ${pct(S.hlMin)} here</span></div>${foot}`;
    return;
  }
  const e = S.active, st = e.st, L = st.L;
  const imgs = [];
  for (let i = 0; i < S.sc[0]; i++) for (let j = 0; j < S.sc[1]; j++) for (let k = 0; k < S.sc[2]; k++) imgs.push([i, j, k]);
  // pairs of every type in every shown cell, starting from atoms inside the cells
  const bonds = [], ends = new Map();
  let truncated = false;
  types.forEach((t) => {
    t.count = 0;
    for (const img of imgs) {
      for (const [i, j, n] of typePairs(e, t.k)) {
        if (bonds.length >= MAX_HIGHLIGHT_PAIRS) { truncated = true; return; }
        const nj = vadd(img, n);
        const a = vadd(st.atoms[i].xyz, imgVec(img, L)), b = vadd(st.atoms[j].xyz, imgVec(nj, L));
        bonds.push({ a, b, color: t.color });
        ends.set(`${i}:${img}`, { s: i, x: a });
        ends.set(`${j}:${nj}`, { s: j, x: b });
        t.count++;
      }
    }
  });
  hlGroup = new THREE.Group();
  const glow = (mat, v) => { mat.emissive = new THREE.Color(0xffffff); mat.emissiveIntensity = v; return mat; };
  const bondMesh = new THREE.InstancedMesh(cylGeo, glow(new THREE.MeshStandardMaterial({ roughness: 0.4 }), 0.12), Math.max(bonds.length, 1));
  bonds.forEach((b, k) => {
    placeCylinder(_dummy, b.a, b.b, S.hlRod);
    _dummy.updateMatrix();
    bondMesh.setMatrixAt(k, _dummy.matrix);
    bondMesh.setColorAt(k, _c.set(b.color));
  });
  bondMesh.count = bonds.length;
  const atoms = [...ends.values()];
  const atomMesh = new THREE.InstancedMesh(sphereGeo, glow(new THREE.MeshStandardMaterial({ roughness: 0.4 }), 0.08), Math.max(atoms.length, 1));
  atoms.forEach((a, k) => {
    const r = atomRadius(st.atoms[a.s]) * 1.2;
    _dummy.position.set(...a.x); _dummy.quaternion.identity(); _dummy.scale.set(r, r, r);
    _dummy.updateMatrix();
    atomMesh.setMatrixAt(k, _dummy.matrix);
    atomMesh.setColorAt(k, _c.set(st.atoms[a.s].color));
  });
  atomMesh.count = atoms.length;
  [bondMesh, atomMesh].forEach((m) => { m.computeBoundingSphere(); hlGroup.add(m); });
  root.add(hlGroup);
  $('hlLegend').innerHTML = head
    + types.map((t) => `<div><i style="background:${t.color}"></i>${esc(t.pair)} <b>${t.d.toFixed(3)} Å</b>
      ${pct(t.share)}${t.negative ? ' <span class="neg">negative</span>' : ''} <span>· ${t.n} pairs per cell</span></div>`).join('')
    + (truncated ? `<div><span>only the first ${MAX_HIGHLIGHT_PAIRS} pairs are drawn</span></div>` : '') + foot;
}
const _dummy = new THREE.Object3D();
function placeCylinder(mesh, a, b, r) {
  _d.set(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
  const len = _d.length() || 1e-6;
  mesh.quaternion.setFromUnitVectors(_Y, _d.multiplyScalar(1 / len));
  mesh.position.set((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2);
  mesh.scale.set(r, len, r);
}

// hover: element, label and occupancy of the atom under the mouse
const ray = new THREE.Raycaster(), mouse = new THREE.Vector2();
renderer.domElement.addEventListener('pointermove', (ev) => {
  const tip = $('tooltip');
  if (ev.buttons || !S.active) { tip.hidden = true; return; }
  const b = renderer.domElement.getBoundingClientRect();
  mouse.set(((ev.clientX - b.left) / b.width) * 2 - 1, -((ev.clientY - b.top) / b.height) * 2 + 1);
  ray.setFromCamera(mouse, camera);
  const hit = ray.intersectObjects(meshes.filter((m) => m.userData.atoms))[0];
  if (!hit) { tip.hidden = true; return; }
  const a = S.active.st.atoms[hit.object.userData.items[hit.instanceId].s];
  tip.innerHTML = `<b>${esc(a.label)}</b> ${esc(a.species)}`;
  tip.hidden = false;
  tip.style.left = `${ev.clientX - b.left + 14}px`; tip.style.top = `${ev.clientY - b.top + 10}px`;
});
renderer.domElement.addEventListener('pointerleave', () => { $('tooltip').hidden = true; });

function loop() {
  if (needsRender) {
    needsRender = false;
    renderer.render(scene, camera);
    if (S.active) renderCompass();
  }
  requestAnimationFrame(loop);
}
requestAnimationFrame(loop);

// ------------------------------------------------------------------ examples from the COD
const CATEGORIES = { inorganic: 'Inorganic', organic: 'Organic', hybrid: 'Hybrid (metal–organic)' };
(async () => {
  try {
    const catalog = await (await fetch('examples/catalog.json')).json();
    const sel = $('examples');
    Object.entries(CATEGORIES).forEach(([key, title]) => {
      const group = document.createElement('optgroup');
      group.label = title;
      catalog.filter((x) => x.category === key).forEach((x) => {
        const o = new Option(`${x.name} – ${x.formula}`, x.file);
        o.title = `COD ${x.cod}: ${x.reference}${x.temperature ? `, ${x.temperature} K` : ''}`;
        group.append(o);
      });
      sel.append(group);
    });
    sel.onchange = async () => {
      const x = catalog.find((c) => c.file === sel.value);
      sel.value = '';
      if (!x) return;
      try {
        const text = await (await fetch(`examples/${x.file}`)).text();
        const e = await loadText(`COD ${x.cod}.cif`, text, x.name);
        if (e) { e.source = x; renderStructInfo(); }
      } catch (err) { fail(err); }
    };
  } catch (err) { console.error(err); }
})();

// ------------------------------------------------------------------ files: button and drag & drop
$('file').onchange = async () => { const files = [...$('file').files]; $('file').value = ''; await openFiles(files); };
let dragDepth = 0;
const hasFiles = (ev) => [...(ev.dataTransfer?.types || [])].includes('Files');
window.addEventListener('dragenter', (ev) => { if (hasFiles(ev)) { dragDepth++; $('dropHint').hidden = false; } });
window.addEventListener('dragleave', (ev) => { if (hasFiles(ev) && --dragDepth <= 0) { dragDepth = 0; $('dropHint').hidden = true; } });
window.addEventListener('dragover', (ev) => { if (hasFiles(ev)) ev.preventDefault(); });
window.addEventListener('drop', (ev) => {
  if (!hasFiles(ev)) return;
  ev.preventDefault();
  dragDepth = 0; $('dropHint').hidden = true;
  openFiles([...ev.dataTransfer.files]);
});

// ------------------------------------------------------------------ info icons
// a small "i" next to headings/controls; its data-tip text appears in a floating box on hover
(() => {
  const tip = document.createElement('div');
  tip.id = 'infotip'; tip.hidden = true; document.body.appendChild(tip);
  const render = (text) => {
    tip.replaceChildren(...text.split('\n').map((line) => {
      const p = document.createElement('div');
      const m = line.match(/^(Ticked|Unticked):\s*/);
      if (m) {
        const b = document.createElement('b');
        b.className = 'tick'; b.textContent = `${m[1]}:`;
        p.append(b, ` ${line.slice(m[0].length)}`);
      } else p.textContent = line;
      return p;
    }));
  };
  const show = (el) => {
    render(el.dataset.tip); tip.hidden = false;
    const r = el.getBoundingClientRect(), w = tip.offsetWidth, h = tip.offsetHeight;
    const x = Math.min(window.innerWidth - w - 8, Math.max(8, r.left + r.width / 2 - w / 2));
    let y = r.bottom + 6;
    if (y + h > window.innerHeight - 8) y = r.top - h - 6;
    tip.style.left = `${x}px`; tip.style.top = `${y}px`;
  };
  document.querySelectorAll('i.info').forEach((el) => { el.tabIndex = 0; });
  document.addEventListener('mouseover', (e) => { const el = e.target.closest?.('i.info'); if (el) show(el); });
  document.addEventListener('mouseout', (e) => { if (e.target.closest?.('i.info')) tip.hidden = true; });
  document.addEventListener('focusin', (e) => { if (e.target.matches?.('i.info')) show(e.target); });
  document.addEventListener('focusout', (e) => { if (e.target.matches?.('i.info')) tip.hidden = true; });
  // icons sit inside <label>s: a click must not toggle the associated checkbox
  document.addEventListener('click', (e) => { if (e.target.closest?.('i.info')) e.preventDefault(); }, true);
})();

// ------------------------------------------------------------------ control states
// controls that cannot have any effect in the current state are greyed out and locked; a label
// is dimmed only when every control in it is locked
function setEnabled(el, on) {
  if (typeof el === 'string') el = $(el);
  if (!el) return;
  el.disabled = !on;
  const lab = el.closest('label');
  if (lab) lab.classList.toggle('off', [...lab.querySelectorAll('input, select, button')].every((c) => c.disabled));
}
function updateControlStates() {
  const e = S.active, has = !!(e && current);
  // viewer
  ['scA', 'scB', 'scC', 'scCube', 'scRmax', 'showBonds', 'boundaryAtoms', 'showCell', 'savePng'].forEach((id) => setEnabled(id, !!e));
  setEnabled('showH', !!e && e.st.atoms.some((a) => a.el === 'H' || a.el === 'D'));
  document.querySelectorAll('[data-view]').forEach((b) => setEnabled(b, !!e));
  // plot
  ['exportCsv', 'exportSvg', 'exportPng', 'modeZoom', 'modePan', 'rFrom', 'rTo'].forEach((id) => setEnabled(id, has));
  setEnabled('resetZoom', has && !!S.xr);
  setEnabled('resetRange', has && !!S.xr);
  // options
  setEnabled('sigma', !S.showHist);
  setEnabled('sigmaRange', !S.showHist);
  const labels = e ? partialLabels() : [], totalOnly = !!CONV[S.conv].totalOnly;
  setEnabled('curvesAll', has && [...(totalOnly ? [] : labels), 'total'].some((l) => S.hidden.has(l)));
  setEnabled('curvesNone', has && !totalOnly && labels.some((l) => !S.hidden.has(l)));
  setEnabled('hlClear', !!S.hl);
}

setBackground();
updatePdf(false);
renderCalcInfo();

// ------------------------------------------------------------------ result of the command line tool
// `xPDFsim file.cif ...` shows its result here: the file and its settings come from the server
const CLI_CONVENTIONS = { G_r: 'G', g_r: 'g', D_r: 'D', GK_r: 'GK', R_r: 'R', rho_r: 'rho', T_r: 'T', N_r: 'N' };
(async () => {
  try {
    const pre = await (await fetch('api/startup')).json();
    if (!pre.text) return;
    const st = pre.settings || {};
    if (st.r_max) $('rMax').value = st.r_max;
    if (st.dr) $('binWidth').value = st.dr;
    if (st.sigma) { $('sigma').value = st.sigma; setSigma(st.sigma, false); }
    setLayer(!!st.histogram);
    if (CLI_CONVENTIONS[st.pdf_type]) S.conv = CLI_CONVENTIONS[st.pdf_type];
    S.rad = st.radiation === 'neutron' ? 'neutron' : 'xray';
    await loadText(pre.file, pre.text);
  } catch (err) { fail(err); }
})();
