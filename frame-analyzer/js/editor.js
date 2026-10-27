/* Custom-model editor: a drawing with mouse tools (select / move, add node,
 * frame, shell, zero-length element, delete) and editable tables of nodes and
 * elements.
 *
 * The custom model (state, saved with the model file):
 *   nodes   [{id, x, y, z, sup: [0|1 per dof] | null, m (t), P (kN, down)}]
 *   frames  [{id, i, j, sec: 'col'|'beam', w (kN/m), ri, rj}]
 *   shells  [{id, n: [4 node ids], t (m), E (MPa), nu, rho (t/m³), q (kN/m²)}]
 *   links   [{id, i, j, dofs: per dof 'free' | 'rigid' | spring name}]
 *   springs [{name, mat: material selector state}]
 *   massFromLoads
 * Coordinates: x, y in 2D (y vertical); x, y, z in 3D (z vertical).
 */
(function (FA) {
  'use strict';
  const PL = FA.plots;
  const TOL = 1e-6;

  const S = {
    tool: 'select', sel: new Set(), pend: [], view: { az: -58, el: 22 }, snap: 0.5, mesh: [2, 2],
    labels: true, showAll: false, msg: '', undo: [], T: null, drag: null,
  };
  let A = null;                    // api: {cm(), dim(), changed(), beamLoad(), springsChanged(), newSpring()}
  let ROOT = null, CANVAS = null, BAR = null, TABLES = null;

  // ------------------------------------------------------------------------
  // Helpers
  // ------------------------------------------------------------------------
  function h(tag, attrs, ...kids) {
    const e = document.createElement(tag);
    for (const k in attrs || {}) {
      const v = attrs[k];
      if (v === undefined || v === null || v === false) continue;
      if (k === 'class') e.className = v;
      else if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
      else e.setAttribute(k, v === true ? '' : v);
    }
    for (const c of kids.flat()) if (c !== null && c !== undefined && c !== false) e.append(c.nodeType ? c : document.createTextNode(String(c)));
    return e;
  }
  const svgEl = PL.svgEl;
  const ndf = dim => (dim === 2 ? 3 : 6);
  const nextId = list => list.reduce((m, o) => Math.max(m, o.id), 0) + 1;
  const r6 = v => +(+v).toFixed(6);

  const SUPPORTS = {
    2: { free: [0, 0, 0], fixed: [1, 1, 1], pinned: [1, 1, 0], roller: [0, 1, 0] },
    3: { free: [0, 0, 0, 0, 0, 0], fixed: [1, 1, 1, 1, 1, 1], pinned: [1, 1, 1, 0, 0, 0], roller: [0, 0, 1, 0, 0, 0] },
  };
  const DOF_NAMES = { 2: ['Ux', 'Uy', 'Rz'], 3: ['Ux', 'Uy', 'Uz', 'Rx', 'Ry', 'Rz'] };
  function supName(dim, sup) {
    if (!sup || !sup.some(Boolean)) return 'free';
    for (const [k, v] of Object.entries(SUPPORTS[dim])) if (v.every((x, i) => !!x === !!sup[i])) return k;
    return 'custom';
  }

  const cm = () => A.cm();
  const dim = () => A.dim();
  const nodeById = id => cm().nodes.find(n => n.id === id);
  const pos = n => (dim() === 2 ? [n.x, n.y] : [n.x, n.y, n.z]);
  const vert = () => (dim() === 2 ? 1 : 2);

  function snapshot() {
    S.undo.push(JSON.stringify(cm()));
    if (S.undo.length > 60) S.undo.shift();
  }
  function undo() {
    if (!S.undo.length) { S.msg = 'Nothing to undo.'; renderBar(); return; }
    const prev = JSON.parse(S.undo.pop()), m = cm();
    for (const k of Object.keys(m)) delete m[k];
    Object.assign(m, prev);
    S.sel.clear(); S.pend = [];
    S.msg = 'Undone.';
    commit(true);
    A.springsChanged();
  }

  /** After a change of the model: rebuild the analysis model and redraw. */
  function commit(tables) {
    A.changed();
    draw();
    renderBar();
    if (tables) renderTables();
  }

  // ------------------------------------------------------------------------
  // Model operations
  // ------------------------------------------------------------------------
  function findNodeAt(xyz) {
    return cm().nodes.find(n => pos(n).every((v, i) => Math.abs(v - xyz[i]) < TOL)) || null;
  }

  function addNode(xyz, merge = true) {
    if (merge) { const f = findNodeAt(xyz); if (f) return f; }
    const n = { id: nextId(cm().nodes), x: r6(xyz[0]), y: r6(xyz[1]), z: dim() === 3 ? r6(xyz[2]) : 0, sup: null, m: 0, P: 0 };
    cm().nodes.push(n);
    return n;
  }

  /** Sine of the member's angle with the horizontal. */
  function slope(a, b) {
    const pa = pos(a), pb = pos(b), v = vert();
    const d = pb.map((x, i) => x - pa[i]), L = Math.hypot(...d);
    return L > 0 ? Math.abs(d[v]) / L : 0;
  }

  /** New frame element: vertical -> column section; horizontal -> beam section with the beam load; else beam, no load. */
  function addFrame(i, j, props) {
    if (i === j) return null;
    const ex = cm().frames.find(f => (f.i === i && f.j === j) || (f.i === j && f.j === i));
    if (ex) { S.msg = `Frame element ${ex.id} already connects these nodes.`; return null; }
    const sn = slope(nodeById(i), nodeById(j));
    const f = Object.assign({ id: nextId(cm().frames), i, j, sec: sn > 0.999 ? 'col' : 'beam', w: sn < 1e-6 ? A.beamLoad() : 0, ri: 0, rj: 0 }, props || {});
    cm().frames.push(f);
    return f;
  }

  /** Split frame f at the given nodes (lying on it); the pieces keep its properties. */
  function splitFrameAt(f, nodes) {
    const pa = pos(nodeById(f.i)), pb = pos(nodeById(f.j));
    const d = pb.map((v, k) => v - pa[k]), L2 = d.reduce((s, v) => s + v * v, 0);
    const along = n => pos(n).reduce((s, v, k) => s + (v - pa[k]) * d[k], 0) / L2;
    const pts = nodes.slice().sort((a, b) => along(a) - along(b)).map(n => n.id);
    const chain = [f.i].concat(pts, [f.j]);
    const list = cm().frames, at = list.indexOf(f);
    const pieces = [];
    for (let k = 0; k < chain.length - 1; k++) {
      pieces.push({ id: k === 0 ? f.id : 0, i: chain[k], j: chain[k + 1], sec: f.sec, w: f.w,
                    ri: k === 0 ? f.ri : 0, rj: k === chain.length - 2 ? f.rj : 0 });
    }
    list.splice(at, 1, pieces[0]);
    for (const p of pieces.slice(1)) { p.id = nextId(list); list.push(p); }
    return pieces;
  }

  /** Nodes lying strictly inside frame f. */
  function nodesOnFrame(f, cands) {
    const pa = pos(nodeById(f.i)), pb = pos(nodeById(f.j));
    const d = pb.map((v, k) => v - pa[k]), L = Math.hypot(...d);
    return cands.filter(n => {
      if (n.id === f.i || n.id === f.j) return false;
      const p = pos(n), r = p.map((v, k) => v - pa[k]);
      const t = r.reduce((s, v, k) => s + v * d[k], 0) / (L * L);
      if (t <= 1e-6 || t >= 1 - 1e-6) return false;
      return Math.hypot(...r.map((v, k) => v - t * d[k])) < 1e-4 * Math.max(1, L);
    });
  }

  function splitFrameEqual(f, k) {
    const pa = pos(nodeById(f.i)), pb = pos(nodeById(f.j));
    const nodes = [];
    for (let s = 1; s < k; s++) nodes.push(addNode(pa.map((v, c) => v + (pb[c] - v) * s / k)));
    return splitFrameAt(f, nodes);
  }

  function polyArea2(P) {             // twice the area, of the projection on the best plane
    if (P[0].length === 2) {
      let a = 0;
      for (let i = 0; i < 4; i++) { const p = P[i], q = P[(i + 1) % 4]; a += p[0] * q[1] - q[0] * p[1]; }
      return Math.abs(a);
    }
    const d13 = P[2].map((v, k) => v - P[0][k]), d24 = P[3].map((v, k) => v - P[1][k]);
    return Math.hypot(...FA.elements.cross(d13, d24));
  }

  function addShellMesh(ids, nx, ny) {
    const P = ids.map(id => pos(nodeById(id)));
    const span = Math.max(...P.map(p => Math.hypot(...p.map((v, k) => v - P[0][k]))));
    // corners must be clicked in order around the panel (no crossing diagonals)
    const cross = [[0, 1, 2, 3], [0, 2, 1, 3], [0, 1, 3, 2]].map(o => polyArea2(o.map(i => P[i])));
    if (cross[0] < 1e-6 * span * span || cross[0] + 1e-9 < Math.max(cross[1], cross[2]) * 0.999) {
      S.msg = 'Click the 4 corners in order around the panel.';
      return 0;
    }
    snapshot();
    const before = new Set(cm().nodes.map(n => n.id));
    const g = [];
    for (let a = 0; a <= nx; a++) {
      g.push([]);
      for (let b = 0; b <= ny; b++) {
        const s = a / nx, t = b / ny;
        const xyz = P[0].map((_, k) => (1 - s) * (1 - t) * P[0][k] + s * (1 - t) * P[1][k] + s * t * P[2][k] + (1 - s) * t * P[3][k]);
        g[a].push((a === 0 && b === 0) ? ids[0] : (a === nx && b === 0) ? ids[1] : (a === nx && b === ny) ? ids[2] :
          (a === 0 && b === ny) ? ids[3] : addNode(xyz).id);
      }
    }
    const last = cm().shells[cm().shells.length - 1];
    const props = last ? { t: last.t, E: last.E, nu: last.nu, rho: last.rho, q: last.q } : { t: 0.2, E: 25000, nu: 0.2, rho: 2.5, q: 0 };
    let n = 0;
    for (let a = 0; a < nx; a++)
      for (let b = 0; b < ny; b++) {
        cm().shells.push(Object.assign({ id: nextId(cm().shells), n: [g[a][b], g[a + 1][b], g[a + 1][b + 1], g[a][b + 1]] }, props));
        n++;
      }
    // frames along the panel edges are split at the new edge nodes
    const fresh = cm().nodes.filter(nd => !before.has(nd.id));
    for (const f of cm().frames.slice()) {
      const on = nodesOnFrame(f, fresh);
      if (on.length) splitFrameAt(f, on);
    }
    return n;
  }

  function ensureSpring() {
    const m = cm();
    if (!m.springs.length) { m.springs.push(A.newSpring('K1')); A.springsChanged(); }
    return m.springs[0].name;
  }

  function addLink(i, j, dofs) {
    const l = { id: nextId(cm().links), i, j, dofs };
    cm().links.push(l);
    return l;
  }

  function defaultDofs(K, bend) {
    const d = dim();
    if (d === 2) return ['rigid', 'rigid', K];
    const out = Array(6).fill('rigid');
    for (const b of (bend || [3, 4])) out[b] = K;
    return out;
  }

  /** zeroLength spring at the end of frame f nearest to node id `end`. */
  function springAtFrameEnd(f, end) {
    snapshot();
    const K = ensureSpring(), n = nodeById(end);
    const nn = addNode(pos(n), false);
    if (f.i === end) f.i = nn.id; else f.j = nn.id;
    let bend = null;
    if (dim() === 3) {                         // global rotation closest to the member's strong axis
      const pa = pos(nodeById(f.i)), pb = pos(nodeById(f.j));
      const d = pb.map((v, k) => v - pa[k]), L = Math.hypot(...d), ex = d.map(v => v / L);
      const vec = Math.abs(ex[2]) > 1 - 1e-9 ? [0, 1, 0] : FA.elements.cross(ex, [0, 0, 1]);
      const ey = FA.elements.unit(FA.elements.cross(vec, ex)), ez = FA.elements.cross(ex, ey);
      let big = 0;
      for (let k = 1; k < 3; k++) if (Math.abs(ez[k]) > Math.abs(ez[big])) big = k;
      bend = [3 + big];
    }
    addLink(end, nn.id, defaultDofs(K, bend));
    S.sel = new Set(['z' + cm().links[cm().links.length - 1].id]);
    S.msg = `Spring inserted at node ${end} of frame element ${f.id} (new node ${nn.id}). Edit its dofs in the zero-length table.`;
  }

  function deleteItems(keys) {
    const m = cm();
    const nodes = new Set(), fr = new Set(), sh = new Set(), lk = new Set();
    for (const k of keys) {
      const id = +k.slice(1);
      ({ n: nodes, f: fr, s: sh, z: lk })[k[0]].add(id);
    }
    m.frames = m.frames.filter(f => !fr.has(f.id) && !nodes.has(f.i) && !nodes.has(f.j));
    m.shells = m.shells.filter(s => !sh.has(s.id) && !s.n.some(n => nodes.has(n)));
    m.links = m.links.filter(l => !lk.has(l.id) && !nodes.has(l.i) && !nodes.has(l.j));
    m.nodes = m.nodes.filter(n => !nodes.has(n.id));
  }

  // ------------------------------------------------------------------------
  // Drawing
  // ------------------------------------------------------------------------
  function projector() {
    if (dim() === 2) return p => [p[0], p[1]];
    const az = S.view.az * Math.PI / 180, el = S.view.el * Math.PI / 180;
    const ex = [-Math.sin(az), Math.cos(az), 0];
    const ey = [-Math.sin(el) * Math.cos(az), -Math.sin(el) * Math.sin(az), Math.cos(el)];
    return p => [p[0] * ex[0] + p[1] * ex[1] + p[2] * ex[2], p[0] * ey[0] + p[1] * ey[1] + p[2] * ey[2]];
  }

  function transform(W, H) {
    const proj = projector(), m = cm();
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const n of m.nodes) {
      const [x, y] = proj(pos(n));
      x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y);
    }
    if (!isFinite(x0)) { x0 = 0; x1 = 10; y0 = 0; y1 = 6; }
    const span = Math.max(x1 - x0, y1 - y0, 1);
    const pad = dim() === 2 ? Math.max(0.12 * span, 2 * S.snap) : 0.06 * span;
    x0 -= pad; x1 += pad; y0 -= pad; y1 += pad;
    const sc = Math.min((W - 20) / (x1 - x0), (H - 20) / (y1 - y0));
    const cx = 10 + (W - 20 - sc * (x1 - x0)) / 2, cy = 10 + (H - 20 - sc * (y1 - y0)) / 2;
    return {
      W, H, sc, x0, x1, y0, y1,
      S: p => { const q = proj(p); return [cx + (q[0] - x0) * sc, cy + (y1 - q[1]) * sc]; },
      inv: (sx, sy) => [(sx - cx) / sc + x0, y1 - (sy - cy) / sc],
    };
  }

  function draw() {
    if (!CANVAS) return;
    const m = cm(), d = dim();
    CANVAS.innerHTML = '';
    const W = Math.max(300, CANVAS.clientWidth || 700), H = Math.round(Math.min(560, Math.max(320, W * 0.58)));
    const T = S.drag && S.drag.T ? S.drag.T : transform(W, H);
    S.T = T;
    const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, width: W, height: H, class: 'plot editor-svg' }, CANVAS);
    const P = new Map(m.nodes.map(n => [n.id, T.S(pos(n))]));
    const pt = id => P.get(id);

    if (d === 2) {                         // snap grid and axes
      const g = svgEl('g', { class: 'ed-grid' }, svg);
      const step = S.snap * T.sc >= 10 ? S.snap : PL.ticks(T.x0, T.x1, Math.round(W / 60)).step;
      if (step > 0) {
        for (let x = Math.ceil(T.x0 / step) * step; x <= T.x1; x += step) {
          const [sx] = T.S([x, 0]);
          svgEl('line', { x1: sx, x2: sx, y1: 0, y2: H, class: Math.abs(x) < 1e-9 ? 'ed-axis' : '' }, g);
        }
        for (let y = Math.ceil(T.y0 / step) * step; y <= T.y1; y += step) {
          const [, sy] = T.S([0, y]);
          svgEl('line', { x1: 0, x2: W, y1: sy, y2: sy, class: Math.abs(y) < 1e-9 ? 'ed-axis' : '' }, g);
        }
      }
      const tk = PL.ticks(T.x0, T.x1, Math.max(2, Math.round(W / 90)));
      for (const v of tk.values) { const [sx] = T.S([v, 0]); svgEl('text', { x: sx + 2, y: H - 4, class: 'tick' }, svg).textContent = +v.toFixed(3); }
      const tv = PL.ticks(T.y0, T.y1, Math.max(2, Math.round(H / 70)));
      for (const v of tv.values) { const [, sy] = T.S([0, v]); svgEl('text', { x: 3, y: sy - 2, class: 'tick' }, svg).textContent = +v.toFixed(3); }
    } else {
      const proj = projector(), g = svgEl('g', { class: 'triad' }, svg), o0 = [36, H - 30];
      for (const [v, lab] of [[[1, 0, 0], 'X'], [[0, 1, 0], 'Y'], [[0, 0, 1], 'Z']]) {
        const q = proj(v), x = o0[0] + q[0] * 22, y = o0[1] - q[1] * 22;
        svgEl('line', { x1: o0[0], y1: o0[1], x2: x, y2: y, class: 'axis' }, g);
        svgEl('text', { x: x + (q[0] >= 0 ? 3 : -9), y: y + (q[1] >= 0 ? -2 : 10), class: 'tick' }, g).textContent = lab;
      }
    }

    for (const s of m.shells) {
      if (!s.n.every(id => P.has(id))) continue;
      svgEl('polygon', { points: s.n.map(id => pt(id).join(',')).join(' '), class: 'ed-shell' + (S.sel.has('s' + s.id) ? ' sel' : '') }, svg);
    }
    for (const f of m.frames) {
      if (!P.has(f.i) || !P.has(f.j)) continue;
      const a = pt(f.i), b = pt(f.j);
      svgEl('line', { x1: a[0], y1: a[1], x2: b[0], y2: b[1], class: `ed-frame ${f.sec}` + (S.sel.has('f' + f.id) ? ' sel' : '') }, svg);
      const L = Math.hypot(b[0] - a[0], b[1] - a[1]) || 1, off = Math.min(9, 0.25 * L);
      for (const [r, p, q] of [[f.ri, a, b], [f.rj, b, a]]) {
        if (!r) continue;
        svgEl('circle', { cx: p[0] + (q[0] - p[0]) / L * off, cy: p[1] + (q[1] - p[1]) / L * off, r: 3.6, class: 'hinge pinned' }, svg);
      }
    }
    // supports
    for (const n of m.nodes) {
      const kind = supName(d, n.sup);
      if (kind === 'free') continue;
      const [x, y] = pt(n.id);
      if (kind === 'fixed') svgEl('rect', { x: x - 6, y: y, width: 12, height: 5, class: 'support' }, svg);
      else if (kind === 'pinned' || kind === 'roller') {
        svgEl('polygon', { points: `${x},${y} ${x - 6},${y + 9} ${x + 6},${y + 9}`, class: 'ed-sup' }, svg);
        if (kind === 'roller') svgEl('line', { x1: x - 7, x2: x + 7, y1: y + 12, y2: y + 12, class: 'ed-sup-line' }, svg);
      } else svgEl('rect', { x: x - 5, y: y + 1, width: 10, height: 8, class: 'ed-sup' }, svg);
    }
    // zero-length elements: marker next to node j
    for (const l of m.links) {
      if (!P.has(l.j)) continue;
      const [x, y] = pt(l.j);
      svgEl('polygon', { points: `${x + 9},${y - 15} ${x + 15},${y - 9} ${x + 9},${y - 3} ${x + 3},${y - 9}`,
        class: 'ed-link' + (S.sel.has('z' + l.id) ? ' sel' : '') }, svg);
      svgEl('line', { x1: x, y1: y, x2: x + 5, y2: y - 5, class: 'ed-link-line' }, svg);
    }
    for (const n of m.nodes) {
      const [x, y] = pt(n.id), on = S.sel.has('n' + n.id), pend = S.pend.includes(n.id);
      svgEl('circle', { cx: x, cy: y, r: on || pend ? 5 : 3.4, class: 'ed-node' + (on ? ' sel' : '') + (pend ? ' pend' : '') }, svg);
    }
    if (S.labels && m.nodes.length + m.frames.length < 600) {
      for (const n of m.nodes) { const [x, y] = pt(n.id); svgEl('text', { x: x + 5, y: y + 13, class: 'ed-label' }, svg).textContent = n.id; }
      for (const f of m.frames) {
        if (!P.has(f.i) || !P.has(f.j)) continue;
        const a = pt(f.i), b = pt(f.j);
        svgEl('text', { x: (a[0] + b[0]) / 2 + 3, y: (a[1] + b[1]) / 2 - 3, class: 'ed-label el' }, svg).textContent = 'F' + f.id;
      }
      for (const s of m.shells) {
        if (!s.n.every(id => P.has(id))) continue;
        const c = [0, 1].map(k => s.n.reduce((t, id) => t + pt(id)[k], 0) / 4);
        svgEl('text', { x: c[0], y: c[1] + 4, class: 'ed-label el', 'text-anchor': 'middle' }, svg).textContent = 'S' + s.id;
      }
    }
    if (d === 3) svgEl('text', { x: W - 8, y: H - 8, class: 'tick hint', 'text-anchor': 'end' }, svg).textContent =
      'drag empty space to rotate · double-click to reset view';
  }

  // ------------------------------------------------------------------------
  // Picking and mouse tools
  // ------------------------------------------------------------------------
  function svgPoint(ev) {
    const svg = CANVAS.querySelector('svg');
    const r = svg.getBoundingClientRect(), k = S.T.W / r.width;
    return [(ev.clientX - r.left) * k, (ev.clientY - r.top) * k];
  }

  function segDist(p, a, b) {
    const dx = b[0] - a[0], dy = b[1] - a[1], L2 = dx * dx + dy * dy || 1;
    const t = Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2));
    return [Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy), t];
  }

  function inPoly(p, pts) {
    let c = false;
    for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
      const [xi, yi] = pts[i], [xj, yj] = pts[j];
      if ((yi > p[1]) !== (yj > p[1]) && p[0] < (xj - xi) * (p[1] - yi) / (yj - yi) + xi) c = !c;
    }
    return c;
  }

  let cycle = { at: null, k: 0 };
  function pick(p) {
    const m = cm(), T = S.T;
    const near = [];
    for (const n of m.nodes) {
      const q = T.S(pos(n)), d = Math.hypot(q[0] - p[0], q[1] - p[1]);
      if (d < 9) near.push([d, n.id]);
    }
    if (near.length) {                       // coincident nodes: repeated clicks cycle through them
      near.sort((a, b) => a[0] - b[0]);
      const same = near.filter(x => x[0] - near[0][0] < 1.5);
      const key = same.map(x => x[1]).join(',');
      cycle = cycle.at === key ? { at: key, k: cycle.k + 1 } : { at: key, k: 0 };
      return { kind: 'n', id: same[cycle.k % same.length][1], many: same.length };
    }
    cycle = { at: null, k: 0 };
    for (const l of m.links) {
      const n = nodeById(l.j);
      if (!n) continue;
      const q = T.S(pos(n));
      if (Math.hypot(q[0] + 9 - p[0], q[1] - 9 - p[1]) < 8) return { kind: 'z', id: l.id };
    }
    let best = null;
    for (const f of m.frames) {
      const a = nodeById(f.i), b = nodeById(f.j);
      if (!a || !b) continue;
      const [dd, t] = segDist(p, T.S(pos(a)), T.S(pos(b)));
      if (dd < 6 && (!best || dd < best.d)) best = { kind: 'f', id: f.id, d: dd, t };
    }
    if (best) return best;
    for (let k = m.shells.length - 1; k >= 0; k--) {
      const s = m.shells[k];
      const pts = s.n.map(id => nodeById(id)).filter(Boolean).map(n => T.S(pos(n)));
      if (pts.length === 4 && inPoly(p, pts)) return { kind: 's', id: s.id };
    }
    return null;
  }

  const snapV = v => r6(Math.round(v / S.snap) * S.snap);
  function worldAt(p) {                      // 2D only
    const [x, y] = S.T.inv(p[0], p[1]);
    return [snapV(x), snapV(y)];
  }

  /** Node under the pointer, or (2D) a new node at the snapped position. */
  function nodeAtClick(hit, p) {
    if (hit && hit.kind === 'n') return hit.id;
    if (dim() !== 2) return null;
    const xyz = worldAt(p), ex = findNodeAt(xyz);
    if (ex) return ex.id;
    return addNode(xyz).id;
  }

  function select(key, add) {
    if (!add) S.sel.clear();
    if (S.sel.has(key) && add) S.sel.delete(key); else S.sel.add(key);
  }

  function capture(ev) { try { CANVAS.setPointerCapture(ev.pointerId); } catch (e) { /* synthetic or ended pointer */ } }

  function onDown(ev) {
    if (ev.button > 0) return;
    const p = svgPoint(ev), hit = pick(p), tool = S.tool;
    S.msg = '';
    if (tool === 'select') {
      if (hit && hit.kind === 'n' && dim() === 2) {
        S.drag = { kind: 'node', id: hit.id, start: p, moved: false, T: S.T, add: ev.shiftKey };
        capture(ev);
      } else if (hit) {
        select(hit.kind + hit.id, ev.shiftKey);
        if (hit.kind === 'n' && hit.many > 1) S.msg = `${hit.many} coincident nodes here: click again to select the next one.`;
        refreshSel();
      } else if (dim() === 3) {
        S.drag = { kind: 'rotate', start: p, az: S.view.az, el: S.view.el, moved: false };
        capture(ev);
      } else { S.drag = { kind: 'clear', start: p }; }
      return;
    }
    if (tool === 'node') {
      if (dim() !== 2) return;
      if (hit && hit.kind === 'n') { select('n' + hit.id); refreshSel(); return; }
      snapshot();
      const n = addNode(worldAt(p));
      S.sel = new Set(['n' + n.id]);
      S.msg = `Node ${n.id} at (${n.x}, ${n.y}).`;
      commit(true);
      return;
    }
    if (tool === 'frame') {
      if (!S.pend.length) snapshot();
      const id = nodeAtClick(hit, p);
      if (id === null) { S.msg = 'Click an existing node (add nodes in the node table in 3D).'; renderBar(); return; }
      if (S.pend.length && S.pend[S.pend.length - 1] === id) { S.pend = []; S.msg = 'Chain finished.'; draw(); renderBar(); return; }
      S.pend.push(id);
      if (S.pend.length === 2) {
        const f = addFrame(S.pend[0], S.pend[1]);
        if (f) { S.msg = `Frame element ${f.id} added (${f.sec === 'col' ? 'column' : 'beam'} section). Click the next node to continue, the same node again or Esc to stop.`; S.pend = [S.pend[1]]; }
        else S.pend = [S.pend[1]];
        commit(true);
        return;
      }
      commit(true);
      return;
    }
    if (tool === 'shell') {
      if (!S.pend.length) snapshot();
      const id = nodeAtClick(hit, p);
      if (id === null) { S.msg = 'Click an existing node.'; renderBar(); return; }
      if (S.pend.includes(id)) { S.msg = 'That corner is already picked.'; renderBar(); return; }
      S.pend.push(id);
      if (S.pend.length === 4) {
        const n = addShellMesh(S.pend, S.mesh[0], S.mesh[1]);
        if (n) S.msg = `${n} shell element${n > 1 ? 's' : ''} added.`;
        S.pend = [];
        commit(true);
      } else { S.msg = `Corner ${S.pend.length} of 4 picked.`; commit(true); }
      return;
    }
    if (tool === 'link') {
      if (hit && hit.kind === 'f' && !S.pend.length) {
        const f = cm().frames.find(x => x.id === hit.id);
        springAtFrameEnd(f, hit.t < 0.5 ? f.i : f.j);
        commit(true);
        return;
      }
      if (!hit || hit.kind !== 'n') { S.msg = 'Click a frame element near one of its ends, or two nodes.'; renderBar(); return; }
      if (!S.pend.length) { S.pend = [hit.id]; S.msg = `Node ${hit.id} picked: click a second node, or the same node again for a spring under its support.`; draw(); renderBar(); return; }
      const a = S.pend[0];
      S.pend = [];
      if (hit.id === a) {
        const n = nodeById(a);
        if (supName(dim(), n.sup) === 'free') { S.msg = 'Node is not supported: click a second node or a frame element end.'; draw(); renderBar(); return; }
        snapshot();
        const K = ensureSpring(), g = addNode(pos(n), false);
        g.sup = n.sup; n.sup = null;
        addLink(g.id, a, defaultDofs(K));
        S.msg = `Support moved to new node ${g.id}; zero-length element added between ${g.id} and ${a}.`;
      } else {
        snapshot();
        addLink(a, hit.id, defaultDofs(ensureSpring()));
        S.msg = `Zero-length element added between nodes ${a} and ${hit.id}.`;
      }
      S.sel = new Set(['z' + cm().links[cm().links.length - 1].id]);
      commit(true);
      return;
    }
    if (tool === 'delete') {
      if (!hit) return;
      snapshot();
      deleteItems([hit.kind + hit.id]);
      S.sel.delete(hit.kind + hit.id);
      S.msg = `${{ n: 'Node', f: 'Frame element', s: 'Shell element', z: 'Zero-length element' }[hit.kind]} ${hit.id} deleted.`;
      commit(true);
    }
  }

  function onMove(ev) {
    const g = S.drag;
    if (!g || g.kind === 'clear') return;
    const p = svgPoint(ev);
    if (Math.hypot(p[0] - g.start[0], p[1] - g.start[1]) > 3 && !g.moved) {
      g.moved = true;
      if (g.kind === 'node') snapshot();
    }
    if (!g.moved) return;
    if (g.kind === 'rotate') {
      S.view.az = g.az - (p[0] - g.start[0]) * 0.5;
      S.view.el = Math.max(-89, Math.min(89, g.el + (p[1] - g.start[1]) * 0.5));
      draw();
    } else if (g.kind === 'node') {
      const n = nodeById(g.id), [x, y] = worldAt(p);
      if (n.x !== x || n.y !== y) { n.x = x; n.y = y; S.msg = `Node ${n.id} → (${x}, ${y})`; draw(); renderBar(); }
    }
  }

  function onUp(ev) {
    const g = S.drag;
    S.drag = null;
    if (!g) return;
    if (g.kind === 'node') {
      if (g.moved) { S.sel = new Set(['n' + g.id]); commit(true); }
      else { select('n' + g.id, g.add); refreshSel(); }
    } else if (g.kind === 'clear' || (g.kind === 'rotate' && !g.moved)) {
      if (S.sel.size) { S.sel.clear(); refreshSel(); }
    } else draw();
  }

  function refreshSel() { draw(); renderBar(); renderTables(); }

  function onKey(ev) {
    if (!ROOT || !ROOT.isConnected || ROOT.offsetParent === null) return;
    const tag = (ev.target.tagName || '').toLowerCase();
    if (tag === 'input' || tag === 'select' || tag === 'textarea') return;
    if (ev.key === 'Escape') { S.pend = []; S.msg = ''; draw(); renderBar(); }
    else if ((ev.key === 'Delete' || ev.key === 'Backspace') && S.sel.size) { ev.preventDefault(); delSelected(); }
    else if ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === 'z') { ev.preventDefault(); undo(); }
  }

  function delSelected() {
    snapshot();
    deleteItems([...S.sel]);
    S.msg = `${S.sel.size} item${S.sel.size > 1 ? 's' : ''} deleted.`;
    S.sel.clear();
    commit(true);
  }

  // ------------------------------------------------------------------------
  // Toolbar, selection bar and tables
  // ------------------------------------------------------------------------
  const TOOLS = [['select', 'Select / move'], ['node', '+ Node'], ['frame', '+ Frame'], ['shell', '+ Shell'], ['link', '+ Zero-length'], ['delete', 'Delete']];
  const HINTS = {
    select: { 2: 'Click to select (Shift adds), drag a node to move it (snaps to the grid).', 3: 'Click to select (Shift adds). Drag empty space to rotate.' },
    node: { 2: 'Click to add a node at the snapped grid point.', 3: '' },
    frame: { 2: 'Click two nodes (or empty grid points) to add a frame element; keep clicking to draw a chain. Esc stops.', 3: 'Click two nodes to add a frame element; keep clicking to chain. Esc stops.' },
    shell: { 2: 'Click the 4 corners of a panel in order. It is meshed into nx × ny shells; frame elements on its edges are split at the new nodes.', 3: 'Click the 4 corners of a panel in order. It is meshed into nx × ny shells; frame elements on its edges are split at the new nodes.' },
    link: { 2: 'Click a frame element near an end to insert a spring there; or click two nodes; or a supported node twice for a spring under the support.', 3: 'Click a frame element near an end to insert a spring there; or click two nodes; or a supported node twice for a spring under the support.' },
    delete: { 2: 'Click a node or an element to delete it (deleting a node deletes its elements). Ctrl+Z undoes.', 3: 'Click a node or an element to delete it (deleting a node deletes its elements). Ctrl+Z undoes.' },
  };

  function numIn(obj, key, o = {}) {
    const inp = h('input', { type: 'number', step: o.step || 'any', value: obj[key], class: o.cls || '' });
    let snap = false;
    inp.addEventListener('focus', () => { snap = false; });
    inp.addEventListener('input', () => {
      const v = Number(inp.value);
      const ok = inp.value !== '' && isFinite(v) && (o.min === undefined || v >= o.min) && (!o.int || Number.isInteger(v));
      inp.classList.toggle('invalid', !ok);
      if (!ok || v === obj[key]) return;
      if (!snap) { snapshot(); snap = true; }
      obj[key] = v;
      if (o.after) o.after();
      A.changed();
      draw();
    });
    inp.addEventListener('change', () => { if (inp.classList.contains('invalid')) { inp.value = obj[key]; inp.classList.remove('invalid'); } });
    return inp;
  }

  function selIn(options, value, onSet) {
    const s = h('select', {});
    for (const [v, l] of options) s.append(h('option', { value: String(v), selected: String(value) === String(v) }, l));
    s.addEventListener('change', () => { snapshot(); onSet(s.value); A.changed(); draw(); renderBar(); });
    return s;
  }

  function nodeRef(obj, key) {
    return numIn(obj, key, { int: true, min: 1, cls: 'id-in', after: () => { if (!nodeById(obj[key])) S.msg = `Node ${obj[key]} does not exist.`; renderBar(); } });
  }

  function renderToolbar(box) {
    const d = dim();
    const tools = h('div', { class: 'seg ed-tools', role: 'radiogroup' });
    for (const [k, l] of TOOLS) {
      if (k === 'node' && d === 3) continue;
      tools.append(h('button', { type: 'button', class: S.tool === k ? 'on' : '', onclick: () => { S.tool = k; S.pend = []; S.msg = ''; render(ROOT, A); } }, l));
    }
    const opts = h('div', { class: 'ed-opts' });
    if (d === 2) {
      const sn = h('input', { type: 'number', min: 0.01, step: 0.05, value: S.snap, class: 'short' });
      sn.addEventListener('input', () => { const v = +sn.value; if (v > 0) { S.snap = v; draw(); } });
      opts.append(h('label', {}, 'Snap (m) ', sn));
    }
    if (S.tool === 'shell') {
      const mk = k => { const i = h('input', { type: 'number', min: 1, max: 20, step: 1, value: S.mesh[k], class: 'short' }); i.addEventListener('input', () => { const v = +i.value; if (Number.isInteger(v) && v >= 1 && v <= 20) S.mesh[k] = v; }); return i; };
      opts.append(h('label', {}, 'Mesh ', mk(0), ' × ', mk(1)));
    }
    const lab = h('input', { type: 'checkbox', checked: S.labels });
    lab.addEventListener('change', () => { S.labels = lab.checked; draw(); });
    opts.append(h('label', { class: 'check inline' }, lab, 'Labels'));
    opts.append(h('button', { type: 'button', class: 'btn small ghost', title: 'Undo (Ctrl+Z)', onclick: undo }, 'Undo'));
    if (d === 3) {
      const xyz = { x: 0, y: 0, z: 0 };
      const go = () => {
        snapshot();
        const n = addNode([xyz.x, xyz.y, xyz.z]);
        S.sel = new Set(['n' + n.id]);
        S.msg = `Node ${n.id} at (${n.x}, ${n.y}, ${n.z}).`;
        commit(true);
      };
      const f = k => { const i = h('input', { type: 'number', step: 'any', value: 0, class: 'short' }); i.addEventListener('input', () => { xyz[k] = +i.value || 0; }); return i; };
      opts.append(h('span', { class: 'ed-add' }, 'Node at x ', f('x'), ' y ', f('y'), ' z ', f('z'),
        h('button', { type: 'button', class: 'btn small', onclick: go }, 'Add')));
    }
    box.append(h('div', { class: 'ed-toolbar' }, tools, opts));
  }

  function renderBar() {
    if (!BAR) return;
    BAR.innerHTML = '';
    const d = dim(), keys = [...S.sel];
    const nodes = keys.filter(k => k[0] === 'n').map(k => nodeById(+k.slice(1))).filter(Boolean);
    const frames = keys.filter(k => k[0] === 'f').map(k => cm().frames.find(f => f.id === +k.slice(1))).filter(Boolean);
    const btn = (label, fn, cls) => h('button', { type: 'button', class: 'btn small ' + (cls || 'ghost'), onclick: fn }, label);
    const apply = fn => () => { snapshot(); fn(); commit(true); };
    if (keys.length) {
      const row = h('div', { class: 'ed-bar' });
      row.append(h('strong', {}, keys.length === 1 ? describe(keys[0]) : `${keys.length} selected`));
      if (nodes.length) {
        row.append(h('span', { class: 'ed-sep' }, 'Support:'));
        for (const k of ['fixed', 'pinned', 'roller', 'free'])
          row.append(btn(k[0].toUpperCase() + k.slice(1), apply(() => { for (const n of nodes) n.sup = k === 'free' ? null : SUPPORTS[d][k].slice(); })));
      }
      if (frames.length) {
        row.append(h('span', { class: 'ed-sep' }, 'Section:'));
        row.append(btn('Column', apply(() => { for (const f of frames) f.sec = 'col'; })));
        row.append(btn('Beam', apply(() => { for (const f of frames) f.sec = 'beam'; })));
        row.append(h('span', { class: 'ed-sep' }, 'Ends:'));
        row.append(btn('Pin both', apply(() => { for (const f of frames) { f.ri = 1; f.rj = 1; } })));
        row.append(btn('Fix both', apply(() => { for (const f of frames) { f.ri = 0; f.rj = 0; } })));
        row.append(btn('Split ÷2', apply(() => { for (const f of frames) splitFrameEqual(f, 2); S.sel.clear(); })));
      }
      row.append(btn('Delete', delSelected, 'danger'));
      row.append(btn('Clear', () => { S.sel.clear(); refreshSel(); }));
      BAR.append(row);
    }
    BAR.append(h('p', { class: 'note ed-hint' }, S.msg || HINTS[S.tool][d] || ''));
  }

  function describe(key) {
    const id = +key.slice(1), d = dim();
    if (key[0] === 'n') { const n = nodeById(id); return n ? `Node ${id} (${pos(n).join(', ')}) · ${supName(d, n.sup)}` : `Node ${id}`; }
    if (key[0] === 'f') {
      const f = cm().frames.find(x => x.id === id);
      if (!f) return `Frame ${id}`;
      const a = nodeById(f.i), b = nodeById(f.j);
      const L = a && b ? Math.hypot(...pos(b).map((v, k) => v - pos(a)[k])) : NaN;
      return `Frame F${id}: ${f.i} → ${f.j}, L = ${L.toFixed(3)} m, ${f.sec === 'col' ? 'column' : 'beam'} section`;
    }
    if (key[0] === 's') return `Shell S${id}`;
    return `Zero-length Z${id}`;
  }

  function table(title, heads, rows, extra) {
    const t = h('table', { class: 'data ed-table' });
    t.append(h('thead', {}, h('tr', {}, heads.map(x => h('th', { class: 'l' }, x)))));
    const tb = h('tbody');
    for (const r of rows) tb.append(r);
    t.append(tb);
    const body = h('div', {}, h('div', { class: 'table-wrap tall' }, t), extra || null);
    const d = h('details', { class: 'exp', open: S['open_' + title] !== false });
    d.addEventListener('toggle', () => { S['open_' + title] = d.open; });
    d.append(h('summary', {}, title), h('div', { class: 'exp-body' }, body));
    return d;
  }

  function idCell(key, label) {
    return h('td', { class: 'l id', title: 'Select', onclick: () => { select(key, false); refreshSel(); } }, label);
  }
  function delCell(key) {
    return h('td', {}, h('button', { type: 'button', class: 'del', title: 'Delete', onclick: () => { snapshot(); deleteItems([key]); S.sel.delete(key); commit(true); } }, '×'));
  }

  function renderTables() {
    if (!TABLES) return;
    TABLES.innerHTML = '';
    const m = cm(), d = dim(), names = DOF_NAMES[d];
    const total = m.nodes.length + m.frames.length + m.shells.length + m.links.length;
    const limited = total > 160 && !S.showAll;
    const show = key => !limited || S.sel.has(key);
    if (limited) TABLES.append(h('p', { class: 'note' }, `Large model (${total} objects): the tables list the selected objects only. `,
      h('button', { type: 'button', class: 'btn small ghost', onclick: () => { S.showAll = true; renderTables(); } }, 'Show all rows')));

    // nodes
    const supOpts = [['free', 'Free'], ['fixed', 'Fixed'], ['pinned', 'Pinned'], ['roller', 'Roller'], ['custom', 'Custom…']];
    const nrows = [];
    for (const n of m.nodes) {
      const key = 'n' + n.id;
      if (!show(key)) continue;
      const kind = supName(d, n.sup);
      const supCell = h('td', { class: 'l' });
      supCell.append(selIn(supOpts, kind, v => {
        if (v === 'custom') n.sup = (n.sup || SUPPORTS[d].free).slice().map(Boolean).map(Number);
        else n.sup = v === 'free' ? null : SUPPORTS[d][v].slice();
        renderTables();
      }));
      if (kind === 'custom') {
        const fx = h('input', { type: 'text', value: n.sup.join(''), class: 'short mono', title: `1 = restrained, order ${names.join(' ')}` });
        fx.addEventListener('change', () => {
          const s = fx.value.replace(/[^01]/g, '');
          if (s.length !== ndf(d)) { fx.classList.add('invalid'); return; }
          snapshot(); fx.classList.remove('invalid');
          n.sup = s.split('').map(Number);
          A.changed(); draw();
        });
        supCell.append(fx);
      }
      const cells = [idCell(key, n.id), h('td', {}, numIn(n, 'x')), h('td', {}, numIn(n, 'y'))];
      if (d === 3) cells.push(h('td', {}, numIn(n, 'z')));
      cells.push(supCell, h('td', {}, numIn(n, 'm', { min: 0 })), h('td', {}, numIn(n, 'P')), delCell(key));
      nrows.push(h('tr', { class: S.sel.has(key) ? 'sel' : '' }, cells));
    }
    const addN = h('button', { type: 'button', class: 'btn small ghost', onclick: () => {
      snapshot();
      const last = m.nodes[m.nodes.length - 1];
      const xyz = last ? pos(last).map((v, k) => (k === 0 ? v + 1 : v)) : [0, 0, 0];
      const n = addNode(xyz, false);
      S.sel = new Set(['n' + n.id]);
      commit(true);
    } }, '+ Add node');
    TABLES.append(table(`Nodes (${m.nodes.length})`, ['Node', 'x (m)', d === 2 ? 'y (m) ↑' : 'y (m)', ...(d === 3 ? ['z (m) ↑'] : []), 'Support', 'Mass (t)', 'Load P (kN) ↓', ''],
      nrows, h('div', { class: 'actions' }, addN, h('span', { class: 'note' }, `Mass is added to the nodal mass from the gravity loads; P is a downward nodal load. Custom supports: 1 = restrained, order ${names.join(' ')}.`))));

    // frames
    const relOpts = d === 2 ? [[0, 'Rigid'], [1, 'Pinned']] : [[0, 'Rigid'], [1, 'Pin major'], [2, 'Pin both']];
    const frows = [];
    for (const f of m.frames) {
      const key = 'f' + f.id;
      if (!show(key)) continue;
      const a = nodeById(f.i), b = nodeById(f.j);
      const L = a && b ? Math.hypot(...pos(b).map((v, k) => v - pos(a)[k])) : NaN;
      frows.push(h('tr', { class: S.sel.has(key) ? 'sel' : '' },
        idCell(key, 'F' + f.id), h('td', {}, nodeRef(f, 'i')), h('td', {}, nodeRef(f, 'j')),
        h('td', { class: 'l' }, selIn([['col', 'Column'], ['beam', 'Beam']], f.sec, v => { f.sec = v; })),
        h('td', {}, numIn(f, 'w')),
        h('td', { class: 'l' }, selIn(relOpts, f.ri || 0, v => { f.ri = +v; })),
        h('td', { class: 'l' }, selIn(relOpts, f.rj || 0, v => { f.rj = +v; })),
        h('td', {}, isFinite(L) ? L.toFixed(3) : '–'),
        h('td', {}, h('button', { type: 'button', class: 'btn small ghost', title: 'Split into 2 elements', onclick: () => { snapshot(); splitFrameEqual(f, 2); commit(true); } }, '÷2')),
        delCell(key)));
    }
    const addF = h('button', { type: 'button', class: 'btn small ghost', onclick: () => {
      snapshot();
      const ns = m.nodes;
      if (ns.length < 2) return;
      const f = addFrame(ns[ns.length - 2].id, ns[ns.length - 1].id);
      if (f) S.sel = new Set(['f' + f.id]);
      commit(true);
    } }, '+ Add frame element');
    TABLES.append(table(`Frame elements (${m.frames.length})`, ['Elem', 'Node i', 'Node j', 'Section', 'w (kN/m) ↓', 'End i', 'End j', 'L (m)', '', ''], frows,
      h('div', { class: 'actions' }, addF, h('span', { class: 'note' }, 'w is a gravity load per metre of member (self-weight, floor load). Pinned ends release the bending moment (major axis, or both axes in 3D).'))));

    // shells
    if (m.shells.length) {
      const srows = [];
      for (const s of m.shells) {
        const key = 's' + s.id;
        if (!show(key)) continue;
        const nn = h('input', { type: 'text', value: s.n.join(' '), class: 'mono nodes-in' });
        nn.addEventListener('change', () => {
          const v = nn.value.trim().split(/[\s,]+/).map(Number);
          if (v.length !== 4 || v.some(x => !Number.isInteger(x)) || new Set(v).size !== 4) { nn.classList.add('invalid'); return; }
          snapshot(); nn.classList.remove('invalid'); s.n = v; A.changed(); draw();
        });
        srows.push(h('tr', { class: S.sel.has(key) ? 'sel' : '' }, idCell(key, 'S' + s.id), h('td', {}, nn),
          h('td', {}, numIn(s, 't', { min: 0 })), h('td', {}, numIn(s, 'E', { min: 0 })), h('td', {}, numIn(s, 'nu', { min: 0 })),
          h('td', {}, numIn(s, 'rho', { min: 0 })), h('td', {}, numIn(s, 'q')), delCell(key)));
      }
      const setAll = h('button', { type: 'button', class: 'btn small ghost', onclick: () => {
        const s0 = m.shells.find(s => S.sel.has('s' + s.id)) || m.shells[0];
        snapshot();
        for (const s of m.shells) Object.assign(s, { t: s0.t, E: s0.E, nu: s0.nu, rho: s0.rho, q: s0.q });
        commit(true);
      } }, 'Copy properties to all shells');
      TABLES.append(table(`Shell elements (${m.shells.length})`, ['Elem', 'Nodes', 't (m)', 'E (MPa)', 'ν', 'ρ (t/m³)', 'q (kN/m²) ↓', ''], srows,
        h('div', { class: 'actions' }, setAll, h('span', { class: 'note' }, (d === 2 ? 'Linear elastic plane-stress quads (in-plane walls / infills).' :
          'Linear elastic flat shells (MITC4 bending + membrane): walls and slabs.') + ' ρ gives the self-weight, q an area gravity load.'))));
    }

    // zero-length
    if (m.links.length) {
      const dofOpts = [['free', 'Free'], ['rigid', 'Rigid']].concat(m.springs.map(s => [s.name, s.name]));
      const lrows = [];
      for (const l of m.links) {
        const key = 'z' + l.id;
        if (!show(key)) continue;
        const cells = [idCell(key, 'Z' + l.id), h('td', {}, nodeRef(l, 'i')), h('td', {}, nodeRef(l, 'j'))];
        for (let k = 0; k < ndf(d); k++) cells.push(h('td', { class: 'l' }, selIn(dofOpts, l.dofs[k] || 'free', v => { l.dofs[k] = v; })));
        cells.push(delCell(key));
        lrows.push(h('tr', { class: S.sel.has(key) ? 'sel' : '' }, cells));
      }
      TABLES.append(table(`Zero-length elements (${m.links.length})`, ['Elem', 'Node i', 'Node j', ...names, ''], lrows,
        h('p', { class: 'note' }, 'Per dof: Free (no connection), Rigid (node j follows node i, as equalDOF) or a spring material from the sidebar (section 7). Springs act on the relative displacement / rotation of node j with respect to node i in global axes.')));
    }
  }

  /** Build the editor into root. api: see A above. */
  function render(root, api) {
    A = api;
    ROOT = root;
    root.innerHTML = '';
    root.classList.add('editor');
    renderToolbar(root);
    BAR = h('div', { class: 'ed-barwrap' });
    CANVAS = h('div', { class: 'card ed-canvas chart' });
    CANVAS.addEventListener('pointerdown', onDown);
    CANVAS.addEventListener('pointermove', onMove);
    CANVAS.addEventListener('pointerup', onUp);
    CANVAS.addEventListener('dblclick', () => { if (dim() === 3) { S.view.az = -58; S.view.el = 22; draw(); } });
    CANVAS.addEventListener('contextmenu', ev => { if (S.pend.length) { ev.preventDefault(); S.pend = []; draw(); renderBar(); } });
    TABLES = h('div', { class: 'ed-tables' });
    root.append(BAR, CANVAS, TABLES);
    if (typeof ResizeObserver !== 'undefined') {
      let lastW = 0;
      new ResizeObserver(() => { const w = CANVAS.clientWidth; if (w && w !== lastW) { lastW = w; draw(); } }).observe(CANVAS);
    }
    renderBar();
    renderTables();
    draw();
  }

  /** Clear selection and history (a different model was loaded). */
  function reset() { S.sel.clear(); S.pend = []; S.undo = []; S.msg = ''; S.showAll = false; }

  /**
   * Custom model from an engine description (renumbered 1..N). springs:
   * [{name, mat}] used for the spring indices of the links.
   */
  function fromDescription(desc, springs, massFromLoads) {
    const d = desc.dim, ids = new Map();
    desc.nodes.forEach((n, k) => ids.set(n.id, k + 1));
    const nodes = desc.nodes.map(n => ({ id: ids.get(n.id), x: n.xyz[0], y: n.xyz[1], z: d === 3 ? n.xyz[2] : 0,
      sup: n.sup && n.sup.some(Boolean) ? n.sup.slice() : null, m: massFromLoads ? 0 : r6(n.m || 0), P: n.P || 0 }));
    const frames = desc.frames.map((f, k) => ({ id: k + 1, i: ids.get(f.i), j: ids.get(f.j), sec: f.sec, w: f.w || 0, ri: f.ri || 0, rj: f.rj || 0 }));
    const links = desc.links.map((l, k) => ({ id: k + 1, i: ids.get(l.i), j: ids.get(l.j),
      dofs: l.dofs.map(t => (typeof t === 'number' ? springs[t].name : t)) }));
    return { nodes, frames, shells: [], links, springs: springs.slice(), massFromLoads: !!massFromLoads };
  }

  document.addEventListener('keydown', onKey);
  FA.editor = { render, reset, draw, fromDescription, supName, SUPPORTS, DOF_NAMES, _state: S };
})(globalThis.FA = globalThis.FA || {});
