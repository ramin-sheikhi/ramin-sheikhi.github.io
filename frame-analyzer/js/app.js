/* User interface (port of the Streamlit app.py): model form, tabs, analyses,
 * results, CSV export and model files.
 */
(function (FA) {
  'use strict';
  const E = FA.engine, M = FA.materials, SEC = FA.sections, R = FA.records, PL = FA.plots;
  const G = E.G, MPA = E.MPA;
  const C1 = 'var(--c1)', C2 = 'var(--c2)', C3 = 'var(--c3)', C4 = 'var(--c4)', MUTED = 'var(--muted)';
  const STORE_KEY = 'frame-analyzer-state-v1';

  const $ = sel => document.querySelector(sel);
  const clone = x => JSON.parse(JSON.stringify(x));
  let uidN = 0;
  const uid = () => 'f' + (++uidN);

  // ========================================================================
  // State
  // ========================================================================
  const matState = family => ({ type: M.DEFAULT[family], vals: {} });
  const sectionDefaults = d => ({
    kind: 'RC rectangle (fiber)', b: d[0], h: d[1], E: 25000, cover: 40, d_bar: 20, n_top: 3, n_bot: 3, n_side: 0,
    conc: matState(M.C), K: 1.2, epscu_core: 0.015, crack: 0.5, rebar: matState(M.S),
    d: d[2], bf: d[3], tw: 11, tf: 19, steel: matState(M.S),
  });

  function defaults() {
    return {
      version: 1,
      dim: 2,
      geo: { n_stories: 3, n_bays: 2, bay: 5, n_bays_y: 2, bay_y: 5, h1: 3.5, h: 3.0, rigid: true, pdelta: true,
             layout: 'regular', unequal: false, bays: [], bays_y: [], heights: [] },
      ele: { type: 'forceBeamColumn', integ: 'Lobatto', n_ip: 5 },
      col: sectionDefaults([400, 400, 300, 300]),
      beam: sectionDefaults([300, 500, 400, 200]),
      conn: { type: 'Rigid', by: 'k', k: 50000, gamma: 0.7, mat: matState(M.J) },
      loads: { w_floor: 25, w_roof: 20, mass_mode: 'gravity', m_floor: 25, m_roof: 20 },
      modal: { n_modes: null },
      po: { pattern: 'Triangular (inverted)', dir: 1, drift: 3, steps: 200 },
      th: { src: 'Synthetic (filtered white noise)', fmt: 'PEER AT2', units: 'g', dt_in: 0.01,
            pga: 0.35, dur: 20, fg: 2.5, zg: 0.6, seed: 1, amp: 0.2, per: 0.5, sdur: 10,
            scale: 1, zeta: 0.05, modes: '1, 3', dt: 0.01, extra: 2, dir: 1 },
      custom: { d2: null, d3: null },     // custom models (see editor.js), one per dimension
    };
  }

  /** Fill missing keys of a (possibly older or hand-edited) state from the defaults. */
  function merge(def, v) {
    if (Array.isArray(def) || def === null || typeof def !== 'object') return v === undefined ? def : v;
    const out = {};
    for (const k of Object.keys(def)) out[k] = merge(def[k], v ? v[k] : undefined);
    if (v && typeof v === 'object') for (const k of Object.keys(v)) if (!(k in out)) out[k] = v[k];
    return out;
  }

  function loadState() {
    try {
      const s = localStorage.getItem(STORE_KEY);
      return s ? merge(defaults(), JSON.parse(s)) : null;
    } catch (e) { return null; }
  }
  function saveState() {
    try { localStorage.setItem(STORE_KEY, JSON.stringify(ST)); } catch (e) { /* storage unavailable */ }
  }

  let ST = loadState() || defaults();
  let UPLOAD = null;                      // {name, text} of the uploaded record

  // Material values of a material selector, filled with the family defaults
  function matValues(ms, family) {
    if (!M.CATALOG[ms.type] || !M.CATALOG[ms.type].families.includes(family)) ms.type = M.DEFAULT[family];
    const v = ms.vals[ms.type] || (ms.vals[ms.type] = {});
    for (const p of M.CATALOG[ms.type].params) if (!(p.name in v)) v[p.name] = clone(M.defaultFor(p, family));
    return v;
  }
  const matSpec = (ms, family) => ({ type: ms.type, v: clone(matValues(ms, family)) });

  function secSpec(s) {
    if (s.kind === 'Elastic rectangle') return { kind: s.kind, b: s.b, h: s.h, E: s.E };
    if (s.kind === 'RC rectangle (fiber)') {
      const conc = matSpec(s.conc, M.C), conf = M.isConfinable(conc);
      return { kind: s.kind, b: s.b, h: s.h, cover: s.cover, d_bar: s.d_bar, n_top: s.n_top, n_bot: s.n_bot,
               n_side: s.n_side, conc, K: conf ? s.K : 1.0, epscu_core: conf ? s.epscu_core : 0.0,
               crack: s.crack, steel: matSpec(s.rebar, M.S) };
    }
    return { kind: s.kind, d: s.d, bf: s.bf, tw: s.tw, tf: s.tf, steel: matSpec(s.steel, M.S) };
  }

  const isCustom = () => ST.geo.layout === 'custom';
  const curCustom = () => ST.custom['d' + ST.dim];

  /** Keep the unequal-spacing lists as long as the bay / story counts. */
  function syncLists() {
    const g = ST.geo;
    const fit = (list, n, def) => { list.length = Math.min(list.length, n); while (list.length < n) list.push(def(list.length)); };
    fit(g.bays, g.n_bays, () => g.bay);
    fit(g.bays_y, g.n_bays_y, () => g.bay_y);
    fit(g.heights, g.n_stories, i => (i === 0 ? g.h1 : g.h));
  }
  const sum = a => a.reduce((x, y) => x + y, 0);
  function spans() {             // [bays X, bays Y] in m
    const g = ST.geo;
    if (g.unequal) { syncLists(); return [g.bays.slice(), g.bays_y.slice()]; }
    return [Array(g.n_bays).fill(g.bay), Array(g.n_bays_y).fill(g.bay_y)];
  }

  function EI_L() {
    const [Eb, , Ib] = E.elasticProps(secSpec(ST.beam));
    const bx = spans()[0];
    return Eb * Ib / (sum(bx) / bx.length);
  }
  function springK() {
    if (ST.conn.by === 'k') return ST.conn.k;
    const g = ST.conn.gamma;
    return 3 * EI_L() * g / (1 - g);
  }
  function springSpec() {
    const sp = matSpec(ST.conn.mat, M.J);
    for (const p of M.CATALOG[sp.type].params) if (p.role === 'E') sp.v[p.name] = springK();
    return sp;
  }

  function totalBeamLength() {
    const g = ST.geo, [bx, by] = spans();
    if (ST.dim === 2) return sum(bx);
    return sum(bx) * (g.n_bays_y + 1) + sum(by) * (g.n_bays + 1);
  }
  function masses() {
    const ld = ST.loads;
    if (ld.mass_mode === 'manual') return [ld.m_floor, ld.m_roof];
    const L = totalBeamLength();
    return [ld.w_floor * L / G, ld.w_roof * L / G];
  }

  /** Engine parameters (the P dict of app.py). */
  function params() {
    const g = ST.geo, dim = ST.dim;
    const [mF, mR] = masses();
    const conn = { type: ST.conn.type };
    if (conn.type === 'Semi-rigid') conn.mat = springSpec();
    const elastic = ST.ele.type === 'elasticBeamColumn';
    const p = {
      dim, n_stories: g.n_stories, n_bays: g.n_bays, bay: g.bay,
      n_bays_y: dim === 3 ? g.n_bays_y : 0, bay_y: dim === 3 ? g.bay_y : g.bay,
      h1: g.h1, h: g.h, rigid_diaphragm: g.rigid, pdelta: g.pdelta, ele_type: ST.ele.type,
      n_ip: elastic ? 5 : ST.ele.n_ip, integ: elastic ? 'Lobatto' : ST.ele.integ,
      col_sec: secSpec(ST.col), beam_sec: secSpec(ST.beam), conn,
      w_floor: ST.loads.w_floor, w_roof: ST.loads.w_roof, m_floor: mF, m_roof: mR,
    };
    if (g.unequal) {
      const [bx, by] = spans();
      p.bays = bx; p.heights = g.heights.slice();
      if (dim === 3) p.bays_y = by;
    }
    if (isCustom() && curCustom()) p.custom = customDescription();
    return p;
  }

  /** Engine description of the custom model (see engine.js). */
  function customDescription() {
    const cm = curCustom(), dim = ST.dim;
    const names = cm.springs.map(s => s.name);
    return {
      dim,
      nodes: cm.nodes.map(n => ({ id: n.id, xyz: dim === 2 ? [n.x, n.y] : [n.x, n.y, n.z], sup: n.sup, m: n.m, P: n.P })),
      frames: cm.frames.map(f => Object.assign({}, f)),
      shells: cm.shells.map(s => Object.assign({}, s, { n: s.n.slice() })),
      links: cm.links.map(l => ({ id: l.id, i: l.i, j: l.j,
        dofs: l.dofs.map(t => (t === 'free' || t === 'rigid' ? t : names.indexOf(t))) })),
      springs: cm.springs.map(s => matSpec(s.mat, M.J)),
      massFromLoads: cm.massFromLoads,
    };
  }

  const newSpring = name => ({ name, mat: matState(M.J) });

  /** Custom model generated from the current regular-grid inputs. */
  function customFromRegular() {
    const layout = ST.geo.layout;
    ST.geo.layout = 'regular';
    const p = params();
    ST.geo.layout = layout;
    const springs = [];
    if (ST.conn.type === 'Semi-rigid') {
      const sp = springSpec();
      springs.push({ name: 'K1', mat: { type: sp.type, vals: { [sp.type]: clone(sp.v) } } });
    }
    return FA.editor.fromDescription(E.regularLayout(p), springs, ST.loads.mass_mode !== 'manual');
  }

  // ========================================================================
  // DOM helpers and widgets
  // ========================================================================
  function h(tag, attrs, ...kids) {
    const e = document.createElement(tag);
    for (const k in attrs || {}) {
      const v = attrs[k];
      if (v === undefined || v === null || v === false) continue;
      if (k === 'class') e.className = v;
      else if (k === 'html') e.innerHTML = v;
      else if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
      else if (k === 'style' && typeof v === 'object') Object.assign(e.style, v);
      else e.setAttribute(k, v === true ? '' : v);
    }
    for (const c of kids.flat()) if (c !== null && c !== undefined && c !== false) e.append(c.nodeType ? c : document.createTextNode(String(c)));
    return e;
  }

  const md = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>')
    .replace(/\*(.+?)\*/g, '<em>$1</em>');

  function label(text, id, help) {
    return h('label', { for: id }, text, help ? h('span', { class: 'help', title: help }, '?') : null);
  }

  /** Number input bound to obj[key]. o: {min, max, step, int, help, onChange} */
  function num(text, obj, key, o = {}) {
    const id = uid();
    const inp = h('input', { type: 'number', id, step: o.step || 'any', min: o.min, max: o.max, value: obj[key],
                             inputmode: 'decimal' });
    inp.addEventListener('input', () => {
      const v = inp.value === '' ? NaN : Number(inp.value);
      const ok = isFinite(v) && (o.min === undefined || v >= o.min) && (o.max === undefined || v <= o.max) &&
                 (!o.int || Number.isInteger(v));
      inp.classList.toggle('invalid', !ok);
      if (ok && v !== obj[key]) { obj[key] = v; (o.onChange || modelChanged)(); }
    });
    inp.addEventListener('change', () => {
      if (inp.classList.contains('invalid')) { inp.value = obj[key]; inp.classList.remove('invalid'); }
    });
    return h('div', { class: 'field' }, label(text, id, o.help), inp);
  }

  function txt(text, obj, key, o = {}) {
    const id = uid();
    const inp = h('input', { type: 'text', id, value: obj[key] });
    inp.addEventListener('input', () => { obj[key] = inp.value; (o.onChange || modelChanged)(); });
    return h('div', { class: 'field' }, label(text, id, o.help), inp);
  }

  /** Select bound to obj[key]; options: array of values or [value, label] pairs. */
  function sel(text, options, obj, key, o = {}) {
    const id = uid();
    const s = h('select', { id });
    for (const opt of options) {
      const [v, l] = Array.isArray(opt) ? opt : [opt, opt];
      s.append(h('option', { value: String(v), selected: String(obj[key]) === String(v) }, l));
    }
    s.addEventListener('change', () => {
      const raw = s.value;
      obj[key] = typeof obj[key] === 'number' ? Number(raw) : raw;
      if (o.onChange) o.onChange(); else { renderSidebar(); modelChanged(); }
    });
    return h('div', { class: 'field' }, label(text, id, o.help), s);
  }

  /** Segmented control (radio). */
  function seg(text, options, obj, key, o = {}) {
    const box = h('div', { class: 'seg' + (o.small ? ' small' : ''), role: 'radiogroup' });
    for (const opt of options) {
      const [v, l] = Array.isArray(opt) ? opt : [opt, opt];
      box.append(h('button', {
        type: 'button', class: obj[key] === v ? 'on' : '', role: 'radio', 'aria-checked': String(obj[key] === v),
        onclick: () => {
          if (obj[key] === v) return;
          obj[key] = v;
          if (o.onChange) o.onChange(); else { renderSidebar(); modelChanged(); }
        },
      }, l));
    }
    return h('div', { class: 'field' }, text ? h('div', { class: 'lbl' }, text, o.help ? h('span', { class: 'help', title: o.help }, '?') : null) : null, box);
  }

  function check(text, obj, key, o = {}) {
    const inp = h('input', { type: 'checkbox', checked: !!obj[key] });
    inp.addEventListener('change', () => { obj[key] = inp.checked; (o.onChange || modelChanged)(); });
    return h('label', { class: 'check' }, inp, text);
  }

  function slider(text, min, max, step, obj, key, o = {}) {
    const id = uid();
    const out = h('output', {}, o.fmt ? o.fmt(obj[key]) : obj[key]);
    const inp = h('input', { type: 'range', id, min, max, step, value: obj[key] });
    inp.addEventListener('input', () => {
      obj[key] = Number(inp.value);
      out.textContent = o.fmt ? o.fmt(obj[key]) : obj[key];
      (o.onChange || modelChanged)();
    });
    return h('div', { class: 'field' }, label(text, id, o.help), h('div', { class: 'slider-wrap' }, inp, out));
  }

  const row = (cols, ...kids) => h('div', { class: 'row', style: { '--cols': cols } }, ...kids);
  const note = t => h('p', { class: 'note', html: md(t) });

  const OPEN = { '1 · Geometry': true, '2 · Element': true, '5 · Beam–column connection': true, '6 · Loads & mass': true };
  function exp(title, build, cls) {
    const d = h('details', { class: 'exp' + (cls ? ' ' + cls : ''), open: !!OPEN[title] });
    d.addEventListener('toggle', () => { OPEN[title] = d.open; });
    const body = h('div', { class: 'exp-body' });
    build(body);
    d.append(h('summary', {}, title), body);
    return d;
  }

  function fmt(v, d = 2) {
    if (v === null || v === undefined || !isFinite(v)) return '–';
    return Number(v).toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d });
  }
  const pct = (v, d = 1) => fmt(100 * v, d) + '%';

  function metric(lab, value, sub) {
    return h('div', { class: 'metric' }, h('div', { class: 'm-label' }, lab), h('div', { class: 'm-value' }, value),
      sub ? h('div', { class: 'm-sub' }, sub) : null);
  }

  function banner(kind, html) { return h('div', { class: 'banner ' + kind, html }); }

  function table(cols, rows, o = {}) {
    const t = h('table', { class: 'data' });
    t.append(h('thead', {}, h('tr', {}, cols.map(c => h('th', { class: c.left ? 'l' : '' }, c.label)))));
    const tb = h('tbody');
    for (const r of rows) tb.append(h('tr', {}, cols.map(c => h('td', { class: c.left ? 'l' : '' }, c.fmt ? c.fmt(r[c.key], r) : r[c.key]))));
    t.append(tb);
    return h('div', { class: 'table-wrap' + (o.tall ? ' tall' : '') }, t);
  }

  function download(name, text, mime = 'text/csv') {
    const blob = new Blob([text], { type: mime });
    const a = h('a', { href: URL.createObjectURL(blob), download: name });
    document.body.append(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  }

  function csv(columns) {        // columns: [[name, array], ...]
    const n = Math.max(...columns.map(c => c[1].length));
    const lines = [columns.map(c => c[0]).join(',')];
    for (let i = 0; i < n; i++) lines.push(columns.map(c => (c[1][i] === undefined ? '' : c[1][i])).join(','));
    return lines.join('\n') + '\n';
  }
  const csvButton = (name, build, labelText = 'Download CSV') =>
    h('button', { type: 'button', class: 'btn small', onclick: () => download(name, build()) }, labelText);

  function chart(opts) {
    const el = h('div', { class: 'card' });
    const c = h('div');
    el.append(c);
    requestAnimationFrame(() => PL.lineChart(c, opts));
    return el;
  }

  // ========================================================================
  // Error reporting
  // ========================================================================
  // (keywords found in the engine log / error, likely cause and fix)
  const HINTS = [
    [['Invalid Args', 'Invalid #args', 'invalid type', 'insufficient args', 'envelope needs', 'could not read',
      'unknown material', 'envelope deformations'],
     '**A material definition was rejected.** Check the material parameters — for a *Custom* material, the ' +
     'name and the number/order of arguments (see the OpenSees manual for that material).'],
    [['Invalid section'],
     '**A section definition was rejected.** Check the section dimensions (cover, flange and web thicknesses).'],
    [['NCV must be greater', 'EigenSOE failed', 'eigen analysis'],
     '**Eigen analysis failed.** Usually more modes were requested than the model has dynamic degrees of ' +
     'freedom, or the mass is zero. Reduce *Number of modes* / the Rayleigh mode *j*, and check that the ' +
     'floor masses are > 0.'],
    [['singular', 'zero pivot', 'LinSOE::solve', 'factorization'],
     '**The stiffness matrix is singular — the model is a mechanism or has a zero stiffness.** Check section ' +
     'dimensions/E (no zeros), the connection type (e.g. pinned beams with a very flexible structure) and the ' +
     'joint spring stiffness. In a custom model also check the supports (enough restraints?), nodes that are ' +
     'connected only by pinned ends, and zero-length elements whose dofs are all Free. The node and dof of the zero ' +
     'pivot are given in the engine messages.'],
    [['failed to converge', 'analyze failed', 'Analysis::analyze', 'did not converge', 'compatible element forces'],
     '**Nonlinear solution did not converge.** Try a smaller analysis dt or more pushover steps, a lower scale ' +
     'factor / target drift, or check the run with *elasticBeamColumn* first: if the elastic model runs, the ' +
     'problem is material/section softening (e.g. concrete crushing, very low εcu, b = 0 for steel).'],
  ];
  const NONCONV_HINT = HINTS[4][1];

  function logBox(log) {
    const n = log ? log.split('\n').length : 0;
    return exp(`Engine messages (${n} line${n === 1 ? "" : "s"})`, b => b.append(h('pre', { class: 'code' }, log || 'The engine printed no messages.')), 'log');
  }

  function failureBlock(what, err, log) {
    const box = h('div');
    const msg = err && err.message ? err.message : String(err);
    box.append(banner('err', `<p>${md(what)} failed: ${md(msg)}</p>`));
    const text = ((log || '') + '\n' + msg).toLowerCase();
    for (const [keys, hint] of HINTS) {
      if (keys.some(k => text.includes(k.toLowerCase()))) { box.append(banner('info', `<p>${md(hint)}</p>`)); break; }
    }
    const lb = logBox(log);
    lb.open = true;
    box.append(lb);
    box.append(exp('JavaScript stack trace', b => b.append(h('pre', { class: 'code' }, (err && err.stack) || msg)), 'log'));
    return box;
  }

  // ========================================================================
  // Sidebar
  // ========================================================================
  const UNITS = { f: ['MPa', 'kN·m'], k: ['MPa', 'kN·m/rad'], e: ['', 'rad'], r: ['', ''] };
  const STEP = { f: 1, k: 1000, e: 0.0001, r: 0.01 };

  function envelopeInputs(body, v, family) {
    const cols = family === M.J ? ['M (kN·m)', 'θ (rad)'] : ['σ (MPa)', 'ε'];
    const editor = (labelText, rows) => {
      const wrap = h('div');
      const draw = () => {
        wrap.innerHTML = '';
        wrap.append(note(labelText));
        const t = h('table', { class: 'edit' });
        t.append(h('tr', {}, h('th', {}, cols[0]), h('th', {}, cols[1]), h('th')));
        rows.forEach((r, i) => {
          const cell = j => {
            const inp = h('input', { type: 'number', step: 'any', value: r[j] });
            inp.addEventListener('input', () => {
              const x = Number(inp.value);
              if (inp.value !== '' && isFinite(x)) { r[j] = x; modelChanged(); }
            });
            return h('td', {}, inp);
          };
          t.append(h('tr', {}, cell(0), cell(1), h('td', {}, h('button', {
            type: 'button', class: 'del', title: 'Remove point',
            onclick: () => { rows.splice(i, 1); draw(); modelChanged(); },
          }, '×'))));
        });
        wrap.append(t);
        if (rows.length < 7) wrap.append(h('button', {
          type: 'button', class: 'btn small ghost',
          onclick: () => { const l = rows[rows.length - 1] || [0, 0]; rows.push([l[0], +(l[1] * 1.5 || 0.01).toPrecision(4)]); draw(); modelChanged(); },
        }, '+ Add point'));
      };
      draw();
      return wrap;
    };
    const env = v.env;
    body.append(editor('Positive envelope, 2–7 points', env.pos));
    body.append(check('Negative envelope = positive (mirrored)', { sym: env.neg === null }, 'sym', {
      onChange: () => { env.neg = env.neg === null ? clone(env.pos) : null; renderSidebar(); modelChanged(); },
    }));
    if (env.neg !== null) body.append(editor('Negative envelope (enter magnitudes)', env.neg));
  }

  /** Material selector + parameter form. eWidget renders the initial-stiffness input (joint springs). */
  function materialInputs(body, ms, family, labelText, eWidget) {
    const names = M.namesFor(family);
    matValues(ms, family);
    body.append(sel(labelText, names, ms, 'type'));
    const mat = M.CATALOG[ms.type];
    const v = matValues(ms, family);
    if (mat.note) body.append(note(mat.note));
    const grid = [];
    for (const prm of mat.params) {
      if (prm.role === 'E' && eWidget) body.append(eWidget());
      else if (prm.kind === 'env') envelopeInputs(body, v, family);
      else if (prm.kind === 'text') body.append(txt(prm.label, v, prm.name, { help: prm.help }));
      else grid.push(prm);
    }
    const cells = grid.map(prm => {
      const unit = UNITS[prm.kind][family === M.J ? 1 : 0];
      return num(unit ? `${prm.label} (${unit})` : prm.label, v, prm.name, { step: STEP[prm.kind], help: prm.help || undefined });
    });
    if (cells.length) body.append(row(2, ...cells));
  }

  function sectionInputs(body, s) {
    const kinds = ['RC rectangle (fiber)', 'Steel I-section (fiber)', 'Elastic rectangle'];
    body.append(sel('Section type', kinds, s, 'kind'));
    if (s.kind === 'Elastic rectangle') {
      body.append(row(2, num('b (mm)', s, 'b', { min: 1, step: 10 }), num('h (mm)', s, 'h', { min: 1, step: 10 })));
      body.append(num('E (MPa)', s, 'E', { min: 1, step: 1000 }));
    } else if (s.kind === 'RC rectangle (fiber)') {
      body.append(row(3, num('b (mm)', s, 'b', { min: 1, step: 10 }), num('h (mm)', s, 'h', { min: 1, step: 10 }),
        num('cover (mm)', s, 'cover', { min: 1, step: 5, help: 'Measured from the face to the bar centre.' })));
      body.append(row(4, num('Ø bar (mm)', s, 'd_bar', { min: 1, step: 1 }),
        num('# top', s, 'n_top', { min: 2, int: true, step: 1 }),
        num('# bottom', s, 'n_bot', { min: 2, int: true, step: 1 }),
        num('# per side', s, 'n_side', { min: 0, int: true, step: 1 })));
      body.append(h('div', { class: 'sub-head' }, 'Concrete'));
      materialInputs(body, s.conc, M.C, 'Concrete material');
      if (M.isConfinable({ type: s.conc.type })) {
        body.append(row(2, num('Core confinement factor K', s, 'K', { min: 1, step: 0.05 }),
          num('Core εcu', s, 'epscu_core', { step: 0.001 })));
      } else body.append(note('Core uses the same material as the cover (no automatic confinement).'));
      body.append(num('Cracked-stiffness factor (elastic element only)', s, 'crack', { min: 0.01, max: 1, step: 0.05,
        help: 'Applied to I (and J in 3D) for elasticBeamColumn.' }));
      body.append(h('div', { class: 'sub-head' }, 'Reinforcement'));
      materialInputs(body, s.rebar, M.S, 'Steel material');
    } else {
      body.append(row(4, num('d (mm)', s, 'd', { min: 1, step: 10 }), num('bf (mm)', s, 'bf', { min: 1, step: 10 }),
        num('tw (mm)', s, 'tw', { min: 0.1, step: 1 }), num('tf (mm)', s, 'tf', { min: 0.1, step: 1 })));
      materialInputs(body, s.steel, M.S, 'Steel material');
    }
  }

  const LIVE = [];          // captions refreshed on every input change
  const LIVE_MODEL = [];    // captions refreshed after every model build

  function renderSidebar() {
    const sb = $('#sidebar'), root = $('#model-form');
    const scroll = sb.scrollTop;
    root.innerHTML = '';
    LIVE.length = 0;
    LIVE_MODEL.length = 0;
    const dim = ST.dim, g = ST.geo;
    const layoutChanged = () => {
      if (isCustom() && !curCustom()) ST.custom['d' + ST.dim] = customFromRegular();
      FA.editor.reset();
      clearResults(); renderSidebar(); renderAll(true); modelChanged();
    };
    root.append(seg('Model type', [[2, '2D frame'], [3, '3D frame']], ST, 'dim', { onChange: layoutChanged }));
    const custom = isCustom();

    root.append(exp('1 · Geometry', b => {
      b.append(seg('Layout', [['regular', 'Regular grid'], ['custom', 'Custom model']], g, 'layout', {
        help: 'Custom model: edit nodes, frame elements, supports, shells and zero-length elements in the Model tab. ' +
              'It starts as a copy of the regular grid.',
        onChange: layoutChanged,
      }));
      if (custom) {
        const cm = curCustom(), cap = note('');
        const upd = () => {
          cap.innerHTML = md(`**Custom ${dim}D model**: ${cm.nodes.length} nodes, ${cm.frames.length} frame, ${cm.shells.length} shell and ` +
            `${cm.links.length} zero-length elements. Edit it in the **Model** tab.`);
        };
        upd();
        LIVE.push(upd);
        b.append(cap);
        b.append(h('div', { class: 'actions' }, h('button', {
          type: 'button', class: 'btn small ghost', title: 'Replace the custom model by the regular grid defined by the inputs of the Regular layout',
          onclick: () => {
            if (!confirm('Replace the custom model by a copy of the regular grid? Your edits will be lost.')) return;
            ST.custom['d' + ST.dim] = customFromRegular();
            layoutChanged();
          },
        }, 'Re-create from regular grid')));
      } else {
        const lists = h('div');
        const drawLists = () => {
          lists.innerHTML = '';
          if (!g.unequal) return;
          syncLists();
          const listInputs = (labelText, arr, prefix) => {
            const cells = arr.map((_, i) => num(`${prefix}${i + 1}`, arr, i, { min: 0.1, step: 0.1 }));
            return h('div', {}, h('div', { class: 'sub-head' }, labelText), h('div', { class: 'list-grid' }, cells));
          };
          lists.append(listInputs(dim === 3 ? 'Bay widths in X (m), left to right' : 'Bay widths (m), left to right', g.bays, 'Bay '));
          if (dim === 3) lists.append(listInputs('Bay widths in Y (m)', g.bays_y, 'Bay '));
          lists.append(listInputs('Story heights (m), bottom to top', g.heights, 'Story '));
        };
        const counts = { min: 1, int: true, step: 1, onChange: () => { drawLists(); modelChanged(); } };
        b.append(num('Stories', g, 'n_stories', Object.assign({ max: 30 }, counts)));
        if (!g.unequal) {
          b.append(row(2, num(dim === 3 ? 'Bays in X' : 'Bays', g, 'n_bays', Object.assign({ max: 20 }, counts)),
            num(dim === 3 ? 'Bay width X (m)' : 'Bay width (m)', g, 'bay', { min: 0.1, step: 0.5 })));
          if (dim === 3) b.append(row(2, num('Bays in Y', g, 'n_bays_y', Object.assign({ max: 20 }, counts)),
            num('Bay width Y (m)', g, 'bay_y', { min: 0.1, step: 0.5 })));
          b.append(row(2, num('1st story h (m)', g, 'h1', { min: 0.1, step: 0.1 }), num('Typical h (m)', g, 'h', { min: 0.1, step: 0.1 })));
        } else {
          b.append(row(2, num(dim === 3 ? 'Bays in X' : 'Bays', g, 'n_bays', Object.assign({ max: 20 }, counts)),
            dim === 3 ? num('Bays in Y', g, 'n_bays_y', Object.assign({ max: 20 }, counts)) : h('div')));
        }
        b.append(check('Unequal bay widths / story heights', g, 'unequal', {
          onChange: () => { if (g.unequal) syncLists(); renderSidebar(); modelChanged(); },
        }));
        b.append(lists);
        drawLists();
      }
      b.append(check('Rigid floor diaphragm' + (dim === 3 ? ' (rigidDiaphragm)' : ' (equalDOF)'), g, 'rigid'));
      b.append(check(custom ? 'P-Delta for columns and inclined members' : 'P-Delta for columns', g, 'pdelta'));
      if (custom) b.append(note('In a custom model the floors are the levels above the lowest node that carry mass and hold a ' +
        'horizontal member or slab; the diaphragm ties the free nodes of each floor.'));
    }));

    root.append(exp('2 · Element', b => {
      b.append(sel('Element type', ['forceBeamColumn', 'dispBeamColumn', 'elasticBeamColumn'], ST.ele, 'type'));
      if (ST.ele.type !== 'elasticBeamColumn') {
        b.append(row(2, sel('Integration', ['Lobatto', 'Legendre', 'Radau', 'NewtonCotes'], ST.ele, 'integ', { onChange: modelChanged }),
          num('Integration pts', ST.ele, 'n_ip', { min: 2, max: 10, int: true, step: 1 })));
      } else b.append(note('Elastic element uses gross E·A·I of the section (× cracked-stiffness factor for RC).'));
    }));

    root.append(exp('3 · Column section', b => {
      if (custom) b.append(note('Assigned to the frame elements with section *Column* (by default the vertical ones).'));
      if (dim === 3) b.append(note('Column depth h is along global X, width b along Y.'));
      sectionInputs(b, ST.col);
    }));
    root.append(exp('4 · Beam section', b => {
      if (custom) b.append(note('Assigned to the frame elements with section *Beam* (by default the non-vertical ones). Depth h is in the vertical plane of the member.'));
      else if (dim === 3) b.append(note('Used for beams in both X and Y. Depth h is vertical.'));
      sectionInputs(b, ST.beam);
    }));

    if (custom) root.append(exp('5 · Beam–column connection', b => {
      b.append(note('In a custom model connections are modelled element by element: **pinned ends** in the frame-element ' +
        'table (or *Pin both* on a selection) and **springs** with zero-length elements (+ Zero-length tool).'));
    }));
    else root.append(exp('5 · Beam–column connection', b => {
      const c = ST.conn;
      b.append(seg('Connection type', ['Rigid', 'Pinned', 'Semi-rigid'], c, 'type', {
        help: 'Applied at both ends of every beam. Column-to-column joints and the base stay continuous / fixed.',
      }));
      if (c.type === 'Pinned') b.append(note('Beam ends transfer shear and axial force only (no strong-axis moment).'));
      else if (c.type === 'Semi-rigid') {
        b.append(note('Moment–rotation spring: forces in kN·m, deformations in rad.'));
        const eWidget = () => {
          const w = h('div');
          w.append(seg('Define rotational stiffness by', [['k', 'Stiffness k'], ['gamma', 'Fixity factor γ']], c, 'by', { small: true }));
          if (c.by === 'k') w.append(num('k (kN·m/rad)', c, 'k', { min: 1, step: 5000 }));
          else w.append(slider('Fixity factor γ', 0.05, 0.99, 0.01, c, 'gamma', {
            help: 'γ = 1 / (1 + 3EI/(kL)); 0 = pinned, 1 = rigid.', fmt: x => x.toFixed(2) }));
          return w;
        };
        materialInputs(b, c.mat, M.J, 'Spring material', eWidget);
        const cap = note('');
        b.append(cap);
        const upd = () => {
          try {
            const k0 = M.initialE(springSpec()), eil = EI_L();
            cap.innerHTML = md(`Initial k = ${fmt(k0, 0)} kN·m/rad · fixity factor γ = ${fmt(1 / (1 + 3 * eil / k0), 2)} · k·L/EI = ${fmt(k0 / eil, 1)}`);
          } catch (e) { cap.textContent = ''; }
        };
        upd();
        LIVE.push(upd);
      }
    }));

    if (custom) root.append(exp('6 · Loads & mass', b => {
      const cm = curCustom();
      b.append(note('Gravity loads are defined per element (frame w, shell ρ and q) and per node (P) in the Model tab.'));
      b.append(num('Load w of new beams (kN/m)', ST.loads, 'w_floor', { step: 1, onChange: saveState }));
      b.append(seg('Seismic mass', [[true, 'Gravity loads / g + nodal'], [false, 'Nodal masses only']], cm, 'massFromLoads', { small: true }));
      const cap = note('');
      const upd = () => {
        if (!MODEL || !MODEL.custom) { cap.textContent = ''; return; }
        const W = MODEL.total_mass * G;
        cap.textContent = `Floors found: ${MODEL.ns} · floor masses: ${MODEL.floor_mass.map(m => fmt(m, 1)).join(', ')} t · ` +
          `total weight ${fmt(W, 0)} kN`;
      };
      upd();
      LIVE_MODEL.push(upd);
      b.append(cap);
    }));
    else root.append(exp('6 · Loads & mass', b => {
      const ld = ST.loads;
      b.append(row(2, num('Floor beam load (kN/m)', ld, 'w_floor', { step: 1 }), num('Roof beam load (kN/m)', ld, 'w_roof', { step: 1 })));
      if (dim === 3) b.append(note('Applied to every beam in X and Y.'));
      b.append(seg('Seismic mass', [['gravity', 'From gravity loads (w·L/g)'], ['manual', 'Manual']], ld, 'mass_mode', { small: true }));
      if (ld.mass_mode === 'manual') {
        b.append(row(2, num('Floor mass (t)', ld, 'm_floor', { min: 0, step: 1 }), num('Roof mass (t)', ld, 'm_roof', { min: 0, step: 1 })));
      } else {
        const cap = note('');
        const upd = () => { const [mf, mr] = masses(); cap.textContent = `Floor: ${fmt(mf, 2)} t · Roof: ${fmt(mr, 2)} t`; };
        upd();
        LIVE.push(upd);
        b.append(cap);
      }
    }));

    if (custom) root.append(exp('7 · Zero-length materials', b => {
      const cm = curCustom();
      b.append(note('Spring materials for the zero-length elements: force in kN or moment in kN·m, deformation in m or rad.'));
      cm.springs.forEach((sp, k) => {
        const box = h('div', { class: 'spring-box' });
        const name = h('input', { type: 'text', value: sp.name, class: 'mono' });
        name.addEventListener('change', () => {
          const v = name.value.trim();
          if (!v || v === 'free' || v === 'rigid' || cm.springs.some(o => o !== sp && o.name === v)) { name.value = sp.name; return; }
          for (const l of cm.links) l.dofs = l.dofs.map(t => (t === sp.name ? v : t));
          sp.name = v;
          renderAll(true); modelChanged();
        });
        const used = cm.links.some(l => l.dofs.includes(sp.name));
        box.append(h('div', { class: 'spring-head' }, h('label', {}, 'Name'), name,
          h('button', { type: 'button', class: 'btn small ghost danger', disabled: used, title: used ? 'Used by zero-length elements' : 'Delete',
            onclick: () => { cm.springs.splice(k, 1); renderSidebar(); renderAll(true); modelChanged(); } }, 'Delete')));
        materialInputs(box, sp.mat, M.J, 'Material');
        b.append(box);
      });
      b.append(h('button', { type: 'button', class: 'btn small ghost', onclick: () => {
        let k = cm.springs.length + 1;
        while (cm.springs.some(o => o.name === 'K' + k)) k++;
        cm.springs.push(newSpring('K' + k));
        renderSidebar(); renderAll(true); modelChanged();
      } }, '+ Add spring material'));
    }));
    sb.scrollTop = scroll;
  }

  // ========================================================================
  // Model build (debounced on every change)
  // ========================================================================
  let MODEL = null, MODEL_ERR = null, BUILDING = false, buildTimer = null, buildToken = 0;
  let RES = { dim: ST.dim };
  let RUN = null;               // {kind, cancel}

  function modelChanged() {
    saveState();
    for (const f of LIVE) f();
    clearTimeout(buildTimer);
    BUILDING = true;
    renderModelStatus();
    buildTimer = setTimeout(buildNow, 500);
  }

  async function buildNow() {
    const token = ++buildToken;
    BUILDING = true;
    renderModelStatus();
    let info = null, err = null;
    try {
      info = await E.buildModelAsync(params(), { cancelled: () => token !== buildToken });
      if (!info) return;                 // superseded by a newer change
    } catch (e) { err = e; }
    if (token !== buildToken) return;
    MODEL = info; MODEL_ERR = err; BUILDING = false;
    renderModelTab();
    renderModal();
    for (const f of LIVE_MODEL) f();
    $('#app-title').textContent = `${ST.dim}D Frame Analyzer`;
  }

  async function freshModel() { return E.buildModelAsync(params()); }

  const lite = info => { const o = Object.assign({}, info); delete o.domain; return o; };

  function clearResults() { RES = { dim: ST.dim }; }

  function span(info) {
    let s = 0;
    for (let k = 0; k < info.dim; k++) {
      const v = info.all_nodes.map(n => info.coords[n][k]);
      s = Math.max(s, Math.max(...v) - Math.min(...v));
    }
    return s || 1;
  }
  function displaced(info, disp, sc) {
    const out = {};
    for (const n in disp) out[n] = info.coords[n].map((c, i) => c + sc * disp[n][i]);
    return out;
  }
  function autoScale(info, disp) {
    let amp = 0;
    for (const n in disp) for (const v of disp[n]) amp = Math.max(amp, Math.abs(v));
    return 0.15 * span(info) / (amp || 1);
  }

  // ========================================================================
  // Model tab
  // ========================================================================
  function renderModelStatus() {
    const st = $('#model-status');
    if (!st) return;
    st.innerHTML = '';
    if (BUILDING) st.append(h('span', { class: 'spinner' }), 'Building model and running gravity analysis…');
    else if (MODEL) st.append('Model up to date.');
  }

  function sectionLines(s) {
    if (s.kind === 'Steel I-section (fiber)')
      return [`d=${s.d}, bf=${s.bf}, tw=${s.tw}, tf=${s.tf} mm · ${s.steel.type}`];
    const lines = [`${s.b} × ${s.h} mm`];
    if (s.kind === 'RC rectangle (fiber)') {
      const n = SEC.rcBarLayout(s).length, [As, rho] = SEC.rcSteelSummary(s);
      lines.push(`${n}Ø${s.d_bar}  (top ${s.n_top}, bottom ${s.n_bot}, ${s.n_side}/side)`);
      lines.push(`As = ${fmt(As, 0)} mm²  ρ = ${fmt(100 * rho, 2)}%`);
      lines.push(`${s.conc.type} / ${s.steel.type}`);
    }
    return lines;
  }

  const editorApi = {
    cm: () => curCustom(), dim: () => ST.dim, changed: () => modelChanged(), beamLoad: () => ST.loads.w_floor,
    springsChanged: () => renderSidebar(), newSpring,
  };

  /** Model tab: the editor (custom models, kept across rebuilds) and the model summary. */
  function renderModelTab(force) {
    const panel = $('#tab-model');
    if (!isCustom()) { panel.innerHTML = ''; renderModelInfo(panel); return; }
    let ed = $('#model-editor'), box = $('#model-info');
    if (force || !ed || ed._cm !== curCustom()) {
      panel.innerHTML = '';
      ed = h('div', { id: 'model-editor' });
      ed._cm = curCustom();
      box = h('div', { id: 'model-info' });
      panel.append(ed, box);
      FA.editor.render(ed, editorApi);
    }
    box.innerHTML = '';
    renderModelInfo(box);
  }

  function renderModelInfo(panel) {
    panel.append(h('div', { class: 'status', id: 'model-status' }));
    renderModelStatus();
    if (MODEL_ERR) { panel.append(failureBlock('Model build / gravity analysis', MODEL_ERR, MODEL_ERR.log)); return; }
    if (!MODEL) { panel.append(h('div', { class: 'empty' }, 'Building the model…')); return; }
    const info = MODEL, p = params(), dim = info.dim;
    const left = h('div', { class: 'card' }), right = h('div');
    const fr = h('div');
    left.append(fr);
    let title;
    if (info.custom) title = `Custom ${dim}D model · ${info.ns} floor${info.ns === 1 ? '' : 's'} with mass`;
    else {
      const plan = dim === 3 ? `${p.n_bays}×${p.n_bays_y}-bay` : `${p.n_bays}-bay`;
      title = `${p.n_stories}-story, ${plan} frame · ${info.conn.toLowerCase()} beam–column joints`;
    }
    requestAnimationFrame(() => PL.frame(fr, info, { title, height: 420 }));

    const rows = [];
    for (const [name, s] of [['Column', p.col_sec], ['Beam', p.beam_sec]]) {
      const sp = SEC.sectionProps(s);
      const r = { member: name, section: s.kind, E: sp.E / MPA, A: sp.A * 1e4, I: sp.Iz * 1e8, Iw: sp.Iy * 1e8, J: sp.J * 1e8, bars: '–', As: null, rho: null };
      if (s.kind === 'RC rectangle (fiber)') {
        const [As, rho] = SEC.rcSteelSummary(s);
        Object.assign(r, { bars: `${SEC.rcBarLayout(s).length}Ø${s.d_bar}`, As, rho: 100 * rho });
      }
      rows.push(r);
    }
    const f0 = v => fmt(v, 0);
    const cols = [{ key: 'member', label: 'Member', left: true }, { key: 'section', label: 'Section', left: true },
      { key: 'E', label: 'E (MPa)', fmt: f0 }, { key: 'A', label: 'A (cm²)', fmt: f0 }, { key: 'I', label: 'I (cm⁴)', fmt: f0 }];
    if (dim === 3) cols.push({ key: 'Iw', label: 'I weak (cm⁴)', fmt: f0 }, { key: 'J', label: 'J (cm⁴)', fmt: f0 });
    cols.push({ key: 'bars', label: 'Bars' }, { key: 'As', label: 'As (mm²)', fmt: f0 }, { key: 'rho', label: 'ρ (%)', fmt: v => fmt(v, 2) });
    right.append(table(cols, rows));
    const W = info.total_mass * G;
    right.append(h('div', { class: 'metrics' }, metric('Total seismic weight', `${fmt(W, 0)} kN`),
      metric('Total height', `${fmt(info.levels[info.levels.length - 1] - info.levels[0], 2)} m`)));
    if (info.custom) {
      right.append(note(`Nodes: ${info.all_nodes.length} · Frame elements: ${info.frames.length} (${info.columns.length} column, ` +
        `${info.beams.length} beam section) · Shells: ${info.shells.length} · Zero-length: ${info.links.length} · ` +
        `Equations: ${info.n_eq} · Dynamic dofs: ${info.max_modes} · Gravity analysis converged ✓`));
      if (info.ns) right.append(note(`Floor levels (m): ${info.levels.slice(1).map(z => fmt(z, 2)).join(', ')} · base at ${fmt(info.levels[0], 2)} m`));
      else right.append(banner('warn', '<p>No level above the base carries mass: modal, pushover and time-history analyses need mass. ' +
        'Add beam loads (w), shell loads or nodal masses.</p>'));
    } else {
      right.append(note(`Nodes: ${info.all_nodes.length} · Columns: ${info.columns.length} · Beams: ${info.beams.length} · ` +
        `Joint springs: ${info.springs.length} · Equations: ${info.n_eq} · Gravity analysis converged ✓`));
    }
    right.append(logBox(info.domain.logText()));
    if (info.custom) panel.append(h('h3', {}, 'Model summary'), h('div', { class: 'grid2' }, left, right));
    else panel.append(h('div', { class: 'grid2' }, left, right));

    panel.append(h('h3', {}, 'Cross-sections'));
    const secs = h('div', { class: 'thirds' });
    for (const [name, s] of [['Column', p.col_sec], ['Beam', p.beam_sec]]) {
      const c = h('div', { class: 'card' });
      PL.section(c, s, name, sectionLines(s));
      secs.append(c);
    }
    panel.append(secs);
    if ([p.col_sec.kind, p.beam_sec.kind].includes('RC rectangle (fiber)'))
      panel.append(note('Dark grey = confined core, light grey = cover concrete, orange = bars. Vertical dimension is the section depth h.'));
  }

  // ========================================================================
  // Running analyses
  // ========================================================================
  function progressBox(textLabel) {
    const fill = h('div', { class: 'fill' });
    const t = h('span', {}, textLabel);
    const stop = h('button', { type: 'button', class: 'btn small danger', onclick: () => { if (RUN) RUN.cancel = true; stop.disabled = true; } }, 'Stop');
    const box = h('div', { class: 'progress' }, h('div', { class: 'bar' }, fill), h('div', { class: 'ptext' }, t, stop));
    box.set = (f, s) => { fill.style.width = `${Math.round(100 * Math.min(1, Math.max(0, f)))}%`; if (s) t.textContent = s; };
    return box;
  }

  function setRunning(kind) {
    RUN = kind ? { kind, cancel: false } : null;
    document.querySelectorAll('.run-btn').forEach(b => { b.disabled = !!RUN; });
  }

  const dirInput = (obj, key) => (ST.dim === 3 ? sel('Direction', [[1, 'X'], [2, 'Y']], obj, key, { onChange: saveState }) : null);

  // ---------------- Modal ----------------
  function renderModal() {
    const panel = $('#tab-modal');
    panel.innerHTML = '';
    const maxModes = MODEL ? Math.max(1, MODEL.max_modes) : 200;
    if (ST.modal.n_modes === null || ST.modal.n_modes > maxModes) ST.modal.n_modes = Math.min(ST.dim === 2 ? 3 : 6, maxModes);
    const run = h('button', { type: 'button', class: 'btn primary run-btn', disabled: !!RUN, onclick: runModal }, 'Run modal analysis');
    panel.append(h('div', { class: 'controls' },
      num('Number of modes', ST.modal, 'n_modes', { min: 1, max: Math.min(200, maxModes), int: true, step: 1, onChange: saveState,
        help: `This model has at most ${maxModes} modes.` }), run));
    const out = h('div', { id: 'modal-out' });
    panel.append(out);
    if (RES.modalErr) out.append(failureBlock('Eigen analysis', RES.modalErr.e, RES.modalErr.log));
    if (!RES.modal) {
      if (!RES.modalErr) out.append(h('div', { class: 'empty' }, 'Run the modal analysis to see periods, mass participation and mode shapes.'));
      return;
    }
    const { modes, info } = RES.modal;
    const dirs = 'xy'.slice(0, info.dim === 2 ? 1 : 2).split('');
    if (modes.some(m => m.lambda < 0))
      out.append(banner('warn', '<p><strong>Negative eigenvalue(s):</strong> the tangent stiffness after gravity is not positive ' +
        'definite, so the structure is unstable under its gravity loads (typically P-Delta with very flexible columns, e.g. ' +
        'pinned beams). Periods of those modes are computed from |λ|.</p>'));
    const cols = [{ key: 'mode', label: 'Mode' }, { key: 'T', label: 'T (s)', fmt: v => fmt(v, 4) },
      { key: 'f', label: 'f (Hz)', fmt: v => fmt(v, 3) }, { key: 'omega', label: 'ω (rad/s)', fmt: v => fmt(v, 3) }];
    for (const a of dirs) cols.push({ key: 'gamma_' + a, label: 'Γ' + a, fmt: v => fmt(v, 3) },
      { key: 'ratio_' + a, label: 'Mass ratio ' + a, fmt: v => pct(v) }, { key: 'cum_' + a, label: 'Cumulative ' + a, fmt: v => pct(v) });
    out.append(table(cols, modes, { tall: true }));
    out.append(h('div', { class: 'actions' }, csvButton('modal.csv', () => csv(
      [['Mode', modes.map(m => m.mode)], ['T (s)', modes.map(m => m.T)], ['f (Hz)', modes.map(m => m.f)], ['omega (rad/s)', modes.map(m => m.omega)]]
        .concat(...dirs.map(a => [['Gamma' + a, modes.map(m => m['gamma_' + a])], ['Mass ratio ' + a, modes.map(m => m['ratio_' + a])],
          ['Cumulative ' + a, modes.map(m => m['cum_' + a])]]))))));
    const grid = h('div', { class: 'modes-grid' });
    const view = { az: -58, el: 22 };
    for (const m of modes.slice(0, 12)) {
      const c = h('div', { class: 'card' }), fr = h('div');
      c.append(fr);
      grid.append(c);
      requestAnimationFrame(() => PL.frame(fr, info, {
        coords: displaced(info, m.shape, autoScale(info, m.shape)), title: `Mode ${m.mode}  T = ${fmt(m.T, 3)} s`,
        height: 260, view: info.dim === 3 ? view : undefined,
      }));
    }
    out.append(grid);
    if (modes.length > 12) out.append(note('Mode shapes are drawn for the first 12 modes.'));
  }

  async function runModal() {
    if (RUN) return;
    setRunning('modal');
    delete RES.modalErr;
    const out = $('#modal-out');
    out.innerHTML = '';
    out.append(h('div', { class: 'status' }, h('span', { class: 'spinner' }), 'Building model and computing modes…'));
    let info = null;
    try {
      info = await freshModel();
      await new Promise(r => setTimeout(r, 0));
      RES.modal = { modes: E.runModal(info, ST.modal.n_modes), info: lite(info) };
    } catch (e) {
      delete RES.modal;
      RES.modalErr = { e, log: e.log || (info ? info.domain.logText() : '') };
    }
    setRunning(null);
    renderModal();
  }

  // ---------------- Pushover ----------------
  function renderPo() {
    const panel = $('#tab-po');
    panel.innerHTML = '';
    const po = ST.po;
    const run = h('button', { type: 'button', class: 'btn primary run-btn', disabled: !!RUN, onclick: runPo }, 'Run pushover');
    panel.append(h('div', { class: 'controls' },
      sel('Lateral load pattern', ['Triangular (inverted)', 'Uniform', 'First mode'], po, 'pattern', { onChange: saveState }),
      dirInput(po, 'dir'),
      num('Target roof drift (%)', po, 'drift', { min: 0.001, step: 0.5, onChange: saveState }),
      num('Steps', po, 'steps', { min: 10, max: 5000, int: true, step: 10, onChange: saveState }), run));
    const out = h('div', { id: 'po-out' });
    panel.append(out);
    if (RES.poErr) out.append(failureBlock('Pushover analysis', RES.poErr.e, RES.poErr.log));
    if (!RES.po) {
      if (!RES.poErr) out.append(h('div', { class: 'empty' }, 'Run a pushover analysis to see the capacity curve, drift profile and deformed shape.'));
      return;
    }
    const { r, info, log, dir } = RES.po;
    const ok = r.status.startsWith('Reached');
    out.append(banner(ok ? 'ok' : 'warn', `<p>${md(r.status)}${info.dim === 3 ? ` · direction ${'XY'[dir - 1]}` : ''}</p>`));
    if (!ok && !r.status.includes('by user')) out.append(banner('info', `<p>${md(NONCONV_HINT)}</p>`));
    const lb = logBox(log);
    if (!ok) lb.open = true;
    out.append(lb);
    const Vmax = Math.max(...r.V);
    out.append(h('div', { class: 'metrics' }, metric('Max base shear', `${fmt(Vmax, 1)} kN`),
      metric('V/W', fmt(Vmax / r.W, 3)), metric('Roof disp at end', `${fmt(1000 * r.roof[r.roof.length - 1], 1)} mm`)));
    const dr = r.drifts[r.drifts.length - 1].map(v => 100 * v);
    const stories = dr.map((_, i) => i + 1);
    out.append(h('div', { class: 'halves' },
      chart({ title: 'Pushover curve', xlabel: 'Roof drift (%)', ylabel: 'Base shear (kN)', xzero: true, yzero: true,
              series: [{ x: r.roof.map(u => 100 * u / r.H), y: r.V, color: C1 }] }),
      chart({ title: 'Drift profile at final step', xlabel: 'Interstory drift (%)', ylabel: 'Story', yticks: stories, xzero: true,
              series: [{ x: dr, y: stories, color: C2, type: 'profile' }] })));
    if (RES.po.sc === undefined) RES.po.sc = Math.min(200, Math.max(1, Math.round(autoScale(info, r.final))));
    const fr = h('div');
    const draw = () => PL.frame(fr, info, { coords: displaced(info, r.final, RES.po.sc), color: C2,
      title: `Deformed shape at final step (×${RES.po.sc.toFixed(0)})`, height: 400, view: RES.po.view || (RES.po.view = { az: -58, el: 22 }) });
    out.append(h('div', { class: 'card', style: { marginTop: '14px' } },
      slider('Deformed-shape scale', 1, 200, 1, RES.po, 'sc', { onChange: draw }), fr));
    requestAnimationFrame(draw);
    out.append(h('div', { class: 'actions' }, csvButton('pushover.csv', () => csv([
      ['roof_disp_m', r.roof], ['roof_drift_pct', r.roof.map(u => 100 * u / r.H)], ['base_shear_kN', r.V]]))));
  }

  async function runPo() {
    if (RUN) return;
    setRunning('po');
    delete RES.poErr;
    const po = ST.po, dir = ST.dim === 3 ? po.dir : 1;
    const out = $('#po-out');
    out.innerHTML = '';
    const bar = progressBox('Building model…');
    out.append(bar);
    let info = null;
    try {
      info = await freshModel();
      bar.set(0, 'Running pushover…');
      const r = await E.runPushover(info, po.pattern, po.drift, po.steps, dir, {
        progress: f => bar.set(f, `Running pushover… ${Math.round(100 * f)}%`),
        cancelled: () => RUN && RUN.cancel,
      });
      RES.po = { r, info: lite(info), log: info.domain.logText(), dir };
    } catch (e) {
      delete RES.po;
      RES.poErr = { e, log: e.log || (info ? info.domain.logText() : '') };
    }
    setRunning(null);
    renderPo();
  }

  // ---------------- Time history ----------------
  function currentRecord() {
    const t = ST.th;
    if (t.src === 'Upload record') {
      if (!UPLOAD) return null;
      const [acc, dt] = R.readRecord(UPLOAD.text, t.fmt, t.dt_in);
      const div = { g: 1, 'm/s²': G, 'cm/s²': 100 * G }[t.units];
      return [acc.map(a => a / div), dt];
    }
    if (t.src.startsWith('Synthetic')) return R.syntheticRecord(t.pga, t.dur, 0.01, t.fg, t.zg, t.seed);
    return R.sineRecord(t.amp, t.per, t.sdur, 0.005);
  }

  let recTimer = null;
  function recordChanged() {
    saveState();
    clearTimeout(recTimer);
    recTimer = setTimeout(renderRecordPreview, 250);
  }

  function renderRecordPreview() {
    const box = $('#rec-preview');
    if (!box) return;
    box.innerHTML = '';
    let rec = null;
    try { rec = currentRecord(); } catch (e) {
      box.append(banner('err', `<p>Could not read record: ${md(e.message)}</p>`),
        note('Check the file format selection (PEER AT2 / single column / time, acc).'));
    }
    const btn = $('#th-run');
    if (btn) btn.disabled = !!RUN || !rec;
    if (!rec) {
      if (ST.th.src === 'Upload record' && !UPLOAD) box.append(note('Choose a record file to continue.'));
      return;
    }
    const [acc, dt] = rec;
    let pga = 0;
    for (const a of acc) pga = Math.max(pga, Math.abs(a));
    const t = Array.from(acc, (_, i) => i * dt);
    box.append(chart({ title: `Input record · PGA = ${fmt(pga, 3)} g · dt = ${+dt.toPrecision(6)} s · ${acc.length} pts`,
      xlabel: 'Time (s)', ylabel: 'Acc (g)', height: 190, series: [{ x: t, y: Array.from(acc), color: MUTED, width: 1 }] }));
  }

  function renderTh() {
    const panel = $('#tab-th');
    panel.innerHTML = '';
    const t = ST.th;
    panel.append(h('h3', { style: { marginTop: 0 } }, 'Ground motion'));
    panel.append(seg('Source', ['Upload record', 'Synthetic (filtered white noise)', 'Sine wave'], t, 'src', {
      onChange: () => { saveState(); renderTh(); },
    }));
    const opts = { onChange: recordChanged };
    if (t.src === 'Upload record') {
      const file = h('input', { type: 'file', accept: '.at2,.txt,.dat,.csv,.acc,.AT2' });
      file.addEventListener('change', () => {
        const f = file.files[0];
        if (!f) return;
        const rd = new FileReader();
        rd.onload = () => { UPLOAD = { name: f.name, text: String(rd.result) }; renderTh(); };
        rd.readAsText(f);
      });
      const fileField = h('div', { class: 'field' }, h('div', { class: 'lbl' }, 'Record file'), file,
        UPLOAD ? note(`Loaded: ${UPLOAD.name}`) : null);
      panel.append(h('div', { class: 'controls' }, fileField,
        sel('Format', ['PEER AT2', 'Single column', 'Time, Acc'], t, 'fmt', opts),
        sel('Acceleration units', ['g', 'm/s²', 'cm/s²'], t, 'units', opts),
        num('Record dt (s) (single column only)', t, 'dt_in', Object.assign({ min: 1e-6, step: 0.001 }, opts))));
    } else if (t.src.startsWith('Synthetic')) {
      panel.append(h('div', { class: 'controls' },
        num('PGA (g)', t, 'pga', Object.assign({ min: 0, step: 0.05 }, opts)),
        num('Duration (s)', t, 'dur', Object.assign({ min: 0.1, step: 1 }, opts)),
        num('Ground freq. fg (Hz)', t, 'fg', Object.assign({ min: 0.01, step: 0.1 }, opts)),
        num('Ground damping ζg', t, 'zg', Object.assign({ min: 0, step: 0.05 }, opts)),
        num('Random seed', t, 'seed', Object.assign({ int: true, step: 1 }, opts))));
    } else {
      panel.append(h('div', { class: 'controls' },
        num('Amplitude (g)', t, 'amp', Object.assign({ step: 0.05 }, opts)),
        num('Period (s)', t, 'per', Object.assign({ min: 0.001, step: 0.1 }, opts)),
        num('Duration (s)', t, 'sdur', Object.assign({ min: 0.01, step: 1 }, opts))));
    }
    panel.append(h('div', { id: 'rec-preview' }));

    panel.append(h('h3', {}, 'Analysis settings'));
    const run = h('button', { type: 'button', class: 'btn primary run-btn', id: 'th-run', disabled: true, onclick: runTh }, 'Run time history');
    const save = { onChange: saveState };
    panel.append(h('div', { class: 'controls' },
      num('Scale factor', t, 'scale', Object.assign({ step: 0.1 }, save)),
      num('Damping ratio ζ', t, 'zeta', Object.assign({ min: 0, step: 0.01 }, save)),
      txt('Rayleigh modes i, j', t, 'modes', save),
      num('Analysis dt (s)', t, 'dt', Object.assign({ min: 1e-5, step: 0.001 }, save)),
      num('Free-vibration time (s)', t, 'extra', Object.assign({ min: 0, step: 0.5 }, save)),
      dirInput(t, 'dir')));
    panel.append(h('div', { class: 'actions' }, run));
    const out = h('div', { id: 'th-out' });
    panel.append(out);
    renderRecordPreview();
    renderThResults();
  }

  function renderThResults() {
    const out = $('#th-out');
    if (!out) return;
    out.innerHTML = '';
    if (RES.thWarn) out.append(banner('warn', `<p>${md(RES.thWarn)}</p>`));
    if (RES.thErr) out.append(failureBlock('Time-history analysis', RES.thErr.e, RES.thErr.log));
    if (!RES.th) return;
    const { r, info, log, dir } = RES.th;
    const ok = r.status === 'Completed';
    out.append(banner(ok ? 'ok' : 'warn', `<p>${md(r.status)}${info.dim === 3 ? ` · direction ${'XY'[dir - 1]}` : ''}</p>`));
    if (!ok && !r.status.includes('by user')) out.append(banner('info', `<p>${md(NONCONV_HINT)}</p>`));
    const lb = logBox(log);
    if (!ok) lb.open = true;
    out.append(lb);
    const [aM, bK] = r.rayleigh;
    out.append(note(`Rayleigh: αM = ${aM.toFixed(4)}, βK = ${bK.toFixed(6)} (T_i = ${r.periods[0].toFixed(3)} s, T_j = ${r.periods[1].toFixed(3)} s)`));
    let iMax = 0;
    r.roof.forEach((v, i) => { if (Math.abs(v) > Math.abs(r.roof[iMax])) iMax = i; });
    const Vabs = Math.max(...r.V.map(Math.abs));
    const nst = r.drifts[0].length;
    const dmax = Array.from({ length: nst }, (_, j) => 100 * Math.max(...r.drifts.map(d => Math.abs(d[j]))));
    out.append(h('div', { class: 'metrics' },
      metric('Peak roof disp', `${fmt(1000 * Math.abs(r.roof[iMax]), 1)} mm`, `t = ${fmt(r.t[iMax], 2)} s`),
      metric('Peak base shear', `${fmt(Vabs, 1)} kN`), metric('V/W', fmt(Vabs / r.W, 3)),
      metric('Max interstory drift', `${fmt(Math.max(...dmax), 2)} %`)));
    const stories = dmax.map((_, i) => i + 1);
    const roofMm = r.roof.map(u => 1000 * u);
    out.append(h('div', { class: 'halves' },
      chart({ title: 'Roof displacement (relative)', xlabel: 'Time (s)', ylabel: 'Roof displacement (mm)', series: [{ x: r.t, y: roofMm, color: C1, width: 1.2 }] }),
      chart({ title: 'Base shear', xlabel: 'Time (s)', ylabel: 'Base shear (kN)', series: [{ x: r.t, y: r.V, color: C2, width: 1.2 }] }),
      chart({ title: 'Base shear vs roof displacement', xlabel: 'Roof displacement (mm)', ylabel: 'Base shear (kN)', series: [{ x: roofMm, y: r.V, color: C3, width: 1 }] }),
      chart({ title: 'Peak drift profile', xlabel: 'Peak interstory drift (%)', ylabel: 'Story', yticks: stories, xzero: true, series: [{ x: dmax, y: stories, color: C4, type: 'profile' }] })));
    const cols = [['time_s', r.t], ['roof_disp_m', r.roof], ['base_shear_kN', r.V]];
    for (let i = 0; i < nst; i++) cols.push([`drift_story${i + 1}`, r.drifts.map(d => d[i])]);
    out.append(h('div', { class: 'actions' }, csvButton('time_history.csv', () => csv(cols))));
  }

  async function runTh() {
    if (RUN) return;
    const t = ST.th;
    let rec;
    try { rec = currentRecord(); } catch (e) { return; }
    if (!rec) return;
    delete RES.thErr; delete RES.thWarn;
    let mi = 1, mj = 3;
    const parts = String(t.modes).split(',').map(s => Number(s.trim()));
    if (parts.length === 2 && parts.every(x => Number.isInteger(x) && x > 0)) [mi, mj] = parts;
    else RES.thWarn = 'Could not read the Rayleigh modes; using 1, 3.';
    setRunning('th');
    const out = $('#th-out');
    out.innerHTML = '';
    const bar = progressBox('Building model…');
    out.append(bar);
    const dir = ST.dim === 3 ? t.dir : 1;
    let info = null;
    try {
      info = await freshModel();
      bar.set(0, 'Running time history…');
      const r = await E.runTimeHistory(info, rec[0], rec[1], t.scale, t.zeta, t.dt, t.extra, [mi, mj], dir, {
        progress: f => bar.set(f, `Running time history… ${Math.round(100 * f)}%`),
        cancelled: () => RUN && RUN.cancel,
      });
      RES.th = { r, info: lite(info), log: info.domain.logText(), dir };
    } catch (e) {
      delete RES.th;
      RES.thErr = { e, log: e.log || (info ? info.domain.logText() : '') };
    }
    setRunning(null);
    renderThResults();
    const btn = $('#th-run');
    if (btn) btn.disabled = false;
  }

  // ========================================================================
  // Page wiring
  // ========================================================================
  function renderAll(force) {
    $('#app-title').textContent = `${ST.dim}D Frame Analyzer`;
    renderModelTab(force);
    renderModal();
    renderPo();
    renderTh();
  }

  function activateTab(name) {
    document.querySelectorAll('.tab').forEach(b => b.classList.toggle('active', b.dataset.tab === name));
    document.querySelectorAll('.tab-panel').forEach(p => p.classList.toggle('active', p.id === 'tab-' + name));
    try { localStorage.setItem(STORE_KEY + '-tab', name); } catch (e) { /* ignore */ }
  }

  function init() {
    document.querySelectorAll('.tab').forEach(b => b.addEventListener('click', () => activateTab(b.dataset.tab)));
    try { const t = localStorage.getItem(STORE_KEY + '-tab'); if (t) activateTab(t); } catch (e) { /* ignore */ }
    $('#btn-menu').addEventListener('click', () => $('#sidebar').classList.toggle('open'));
    $('#btn-save').addEventListener('click', () => download('frame-model.json', JSON.stringify(ST, null, 2), 'application/json'));
    $('#btn-open').addEventListener('click', () => $('#file-open').click());
    $('#file-open').addEventListener('change', ev => {
      const f = ev.target.files[0];
      if (!f) return;
      const rd = new FileReader();
      rd.onload = () => {
        try {
          const s = JSON.parse(String(rd.result));
          if (!s || (s.dim !== 2 && s.dim !== 3) || !s.geo) throw new Error('not a model file');
          ST = merge(defaults(), s);
          clearResults();
          FA.editor.reset();
          renderSidebar(); renderAll(true); modelChanged();
        } catch (e) { alert('Could not open the model file: ' + e.message); }
        ev.target.value = '';
      };
      rd.readAsText(f);
    });
    $('#btn-reset').addEventListener('click', () => {
      if (!confirm('Restore the default model? Your current inputs will be replaced.')) return;
      ST = defaults();
      clearResults();
      FA.editor.reset();
      renderSidebar(); renderAll(true); modelChanged();
    });
    renderSidebar();
    renderAll();
    modelChanged();
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();

  FA.app = { get state() { return ST; }, params, customFromRegular, get model() { return MODEL; }, get results() { return RES; } };
})(globalThis.FA = globalThis.FA || {});
