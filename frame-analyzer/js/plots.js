/* SVG drawing: line charts, frame (2D / rotatable 3D) and cross-section
 * drawings. Colours come from CSS variables so light and dark themes work.
 */
(function (FA) {
  'use strict';
  const NS = 'http://www.w3.org/2000/svg';

  function svgEl(tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs || {}) if (attrs[k] !== undefined && attrs[k] !== null) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }

  function text(parent, x, y, str, attrs) {
    const t = svgEl('text', Object.assign({ x, y }, attrs || {}), parent);
    t.textContent = str;
    return t;
  }

  // Re-render registered drawings when their container is resized
  const ro = typeof ResizeObserver !== 'undefined' ? new ResizeObserver(entries => {
    for (const en of entries) {
      const el = en.target;
      const w = Math.round(en.contentRect.width);
      if (el._render && w && w !== el._lastW) { el._lastW = w; el._render(); }
    }
  }) : null;

  function mount(el, render) {
    el._render = render;
    el._lastW = el.clientWidth;
    render();
    if (ro && !el._observed) { ro.observe(el); el._observed = true; }
  }

  // ------------------------------------------------------------------------
  // Number formatting and ticks
  // ------------------------------------------------------------------------
  function niceStep(span, n) {
    const raw = span / Math.max(1, n);
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const r = raw / mag;
    return (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
  }

  function ticks(min, max, n) {
    if (!(max > min)) { max = min + 1; }
    const step = niceStep(max - min, n);
    const out = [];
    for (let v = Math.ceil(min / step - 1e-9) * step; v <= max + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    return { values: out, step };
  }

  function fmtTick(v, step) {
    const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
    if (Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0 && d > 6)) return v.toExponential(1);
    return v.toLocaleString('en-US', { minimumFractionDigits: Math.min(d, 6), maximumFractionDigits: Math.min(d, 6) });
  }

  // ------------------------------------------------------------------------
  // Line chart
  // ------------------------------------------------------------------------
  /**
   * opts: { title, xlabel, ylabel, height, series: [{x, y, color, width, markers,
   *         dash, label}], xmin, xmax, ymin, ymax, yticks (array), xzero, yzero }
   * series.type 'profile': y = story numbers, x = values; drawn as a story profile.
   */
  function lineChart(el, opts) {
    mount(el, () => drawLineChart(el, opts));
  }

  function drawLineChart(el, o) {
    el.innerHTML = '';
    el.classList.add('chart');
    const W = Math.max(280, el.clientWidth || 560), H = o.height || 300;
    const m = { l: 62, r: 14, t: o.title ? 30 : 12, b: 42 };
    const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, width: W, height: H, class: 'plot', role: 'img',
      'aria-label': o.title || 'chart' }, el);
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (const s of o.series) {
      for (let i = 0; i < s.x.length; i++) {
        const x = s.x[i], y = s.y[i];
        if (!isFinite(x) || !isFinite(y)) continue;
        if (x < xmin) xmin = x; if (x > xmax) xmax = x;
        if (y < ymin) ymin = y; if (y > ymax) ymax = y;
      }
      if (s.type === 'profile') { ymin = Math.min(ymin, 0.5); ymax = Math.max(ymax, s.y.length + 0.5); }
    }
    if (!isFinite(xmin)) { xmin = 0; xmax = 1; ymin = 0; ymax = 1; }
    if (o.xzero) xmin = Math.min(0, xmin);
    if (o.yzero) ymin = Math.min(0, ymin);
    if (o.xmin !== undefined) xmin = o.xmin;
    if (o.xmax !== undefined) xmax = o.xmax;
    if (o.ymin !== undefined) ymin = o.ymin;
    if (o.ymax !== undefined) ymax = o.ymax;
    if (xmax - xmin < 1e-12) { xmax += 0.5 + Math.abs(xmax) * 0.05; xmin -= o.xzero && xmin === 0 ? 0 : 0.5 + Math.abs(xmin) * 0.05; }
    if (ymax - ymin < 1e-12) { ymax += 0.5 + Math.abs(ymax) * 0.05; ymin -= o.yzero && ymin === 0 ? 0 : 0.5 + Math.abs(ymin) * 0.05; }
    const padY = (ymax - ymin) * 0.05, padX = (xmax - xmin) * 0.02;
    if (o.ymin === undefined && !(o.yzero && ymin === 0)) ymin -= padY;
    if (o.ymax === undefined) ymax += padY;
    if (o.xmin === undefined && !(o.xzero && xmin === 0)) xmin -= padX;
    if (o.xmax === undefined) xmax += padX;
    const pw = W - m.l - m.r, ph = H - m.t - m.b;
    const X = v => m.l + (v - xmin) / (xmax - xmin) * pw;
    const Y = v => m.t + ph - (v - ymin) / (ymax - ymin) * ph;

    if (o.title) text(svg, m.l, 18, o.title, { class: 'plot-title' });
    const g = svgEl('g', {}, svg);
    const xt = ticks(xmin, xmax, Math.max(3, Math.round(pw / 90)));
    const yt = o.yticks ? { values: o.yticks, step: 1 } : ticks(ymin, ymax, Math.max(3, Math.round(ph / 55)));
    for (const v of xt.values) {
      if (v < xmin - 1e-12 || v > xmax + 1e-12) continue;
      svgEl('line', { x1: X(v), x2: X(v), y1: m.t, y2: m.t + ph, class: 'grid' }, g);
      text(g, X(v), m.t + ph + 16, fmtTick(v, xt.step), { class: 'tick', 'text-anchor': 'middle' });
    }
    for (const v of yt.values) {
      if (v < ymin - 1e-12 || v > ymax + 1e-12) continue;
      svgEl('line', { x1: m.l, x2: m.l + pw, y1: Y(v), y2: Y(v), class: 'grid' }, g);
      text(g, m.l - 6, Y(v) + 4, o.yticks ? String(v) : fmtTick(v, yt.step), { class: 'tick', 'text-anchor': 'end' });
    }
    svgEl('line', { x1: m.l, x2: m.l, y1: m.t, y2: m.t + ph, class: 'axis' }, g);
    svgEl('line', { x1: m.l, x2: m.l + pw, y1: m.t + ph, y2: m.t + ph, class: 'axis' }, g);
    if (o.xlabel) text(g, m.l + pw / 2, H - 6, o.xlabel, { class: 'axis-label', 'text-anchor': 'middle' });
    if (o.ylabel) {
      const t = text(g, 0, 0, o.ylabel, { class: 'axis-label', 'text-anchor': 'middle' });
      t.setAttribute('transform', `translate(14 ${m.t + ph / 2}) rotate(-90)`);
    }
    const clipId = 'c' + Math.random().toString(36).slice(2);
    const cp = svgEl('clipPath', { id: clipId }, svg);
    svgEl('rect', { x: m.l, y: m.t - 2, width: pw + 2, height: ph + 4 }, cp);
    const plot = svgEl('g', { 'clip-path': `url(#${clipId})` }, svg);
    const pts = [];
    o.series.forEach((s, si) => {
      const style = `stroke:${s.color || 'var(--c1)'};stroke-width:${s.width || 2}` + (s.dash ? `;stroke-dasharray:${s.dash}` : '');
      let d = '';
      if (s.type === 'profile') {
        s.x.forEach((v, i) => {
          const y0 = Y(s.y[i] - 0.5), y1 = Y(s.y[i] + 0.5);
          d += (i ? 'L' : 'M') + X(v).toFixed(1) + ' ' + y0.toFixed(1) + 'L' + X(v).toFixed(1) + ' ' + y1.toFixed(1);
        });
      } else {
        const n = s.x.length;
        const stride = n > 4 * pw ? Math.floor(n / (2 * pw)) : 1;     // decimate very long series
        for (let i = 0; i < n; i += stride) d += (i ? 'L' : 'M') + X(s.x[i]).toFixed(1) + ' ' + Y(s.y[i]).toFixed(1);
        if (stride > 1 && (n - 1) % stride) d += 'L' + X(s.x[n - 1]).toFixed(1) + ' ' + Y(s.y[n - 1]).toFixed(1);
      }
      if (s.line !== false) svgEl('path', { d, class: 'series', style }, plot);
      if (s.markers || s.type === 'profile') {
        s.x.forEach((v, i) => svgEl('circle', { cx: X(v), cy: Y(s.y[i]), r: 4, style: `fill:${s.color || 'var(--c1)'}` }, plot));
      }
      for (let i = 0; i < s.x.length; i++) pts.push([s.x[i], s.y[i], si]);
    });
    // hover read-out
    const hover = svgEl('g', { class: 'hover', style: 'display:none' }, svg);
    const dot = svgEl('circle', { r: 5, class: 'hover-dot' }, hover);
    const tip = document.createElement('div');
    tip.className = 'chart-tip';
    el.appendChild(tip);
    svg.addEventListener('mousemove', ev => {
      const r = svg.getBoundingClientRect(), k = W / r.width;
      const px = (ev.clientX - r.left) * k, py = (ev.clientY - r.top) * k;
      if (px < m.l || px > m.l + pw || py < m.t || py > m.t + ph) { hover.style.display = 'none'; tip.style.display = 'none'; return; }
      let best = null, bd = Infinity;
      const step = Math.max(1, Math.floor(pts.length / 4000));
      for (let i = 0; i < pts.length; i += step) {
        const q = pts[i];
        const dx = X(q[0]) - px, dy = Y(q[1]) - py, dd = dx * dx + dy * dy;
        if (dd < bd) { bd = dd; best = q; }
      }
      if (!best || bd > 900) { hover.style.display = 'none'; tip.style.display = 'none'; return; }
      hover.style.display = '';
      dot.setAttribute('cx', X(best[0])); dot.setAttribute('cy', Y(best[1]));
      dot.setAttribute('style', `stroke:${o.series[best[2]].color || 'var(--c1)'}`);
      const f = v => (Math.abs(v) >= 1000 ? v.toLocaleString('en-US', { maximumFractionDigits: 1 }) : Number(v.toPrecision(4)).toString());
      tip.textContent = `${o.xlabel || 'x'}: ${f(best[0])}   ${o.ylabel || 'y'}: ${f(best[1])}`;
      tip.style.display = 'block';
      tip.style.left = Math.min(r.width - 10, (X(best[0]) / k) + 12) + 'px';
      tip.style.top = (Y(best[1]) / k - 34) + 'px';
    });
    svg.addEventListener('mouseleave', () => { hover.style.display = 'none'; tip.style.display = 'none'; });
  }

  // ------------------------------------------------------------------------
  // Frame drawing
  // ------------------------------------------------------------------------
  /**
   * opts: { coords (displaced positions, optional), color, title, height,
   *         view: {az, el} shared object for 3D (mutated by dragging) }
   */
  function frame(el, info, opts = {}) {
    if (info.dim === 3 && !opts.view) opts.view = { az: -58, el: 22 };
    el._frame = { info, opts };
    const redraw = () => drawFrame(el, el._frame.info, el._frame.opts);
    mount(el, redraw);
    if (!el._drag) {               // rotation handlers (act only on 3D drawings)
      el._drag = true;
      let start = null;
      const view = () => (el._frame.info.dim === 3 ? el._frame.opts.view : null);
      const down = (x, y) => {
        const v = view();
        if (!v) return;
        start = { x, y, az: v.az, el: v.el };
        el.classList.add('dragging');
      };
      const move = (x, y) => {
        const v = view();
        if (!start || !v) return;
        v.az = start.az - (x - start.x) * 0.5;
        v.el = Math.max(-89, Math.min(89, start.el + (y - start.y) * 0.5));
        redraw();
      };
      const up = () => { start = null; el.classList.remove('dragging'); };
      el.addEventListener('mousedown', e => { if (view()) { down(e.clientX, e.clientY); e.preventDefault(); } });
      window.addEventListener('mousemove', e => move(e.clientX, e.clientY));
      window.addEventListener('mouseup', up);
      el.addEventListener('touchstart', e => { const t = e.touches[0]; down(t.clientX, t.clientY); }, { passive: true });
      el.addEventListener('touchmove', e => { if (start) { const t = e.touches[0]; move(t.clientX, t.clientY); e.preventDefault(); } }, { passive: false });
      el.addEventListener('touchend', up);
      el.addEventListener('dblclick', () => { const v = view(); if (v) { v.az = -58; v.el = 22; redraw(); } });
    }
  }

  function drawFrame(el, info, o) {
    el.innerHTML = '';
    el.classList.add('chart', 'frame-view');
    const W = Math.max(220, el.clientWidth || 520), H = o.height || Math.round(Math.min(460, Math.max(260, W * 0.72)));
    const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, width: W, height: H, class: 'plot', role: 'img',
      'aria-label': o.title || 'frame drawing' }, el);
    const top = o.title ? 26 : 8;
    if (o.title) text(svg, 10, 17, o.title, { class: 'plot-title small' });
    const d3 = info.dim === 3;
    const base = info.coords;
    let proj;
    if (d3) {
      const az = o.view.az * Math.PI / 180, elv = o.view.el * Math.PI / 180;
      const ex = [-Math.sin(az), Math.cos(az), 0];
      const ey = [-Math.sin(elv) * Math.cos(az), -Math.sin(elv) * Math.sin(az), Math.cos(elv)];
      proj = p => [p[0] * ex[0] + p[1] * ex[1] + p[2] * ex[2], p[0] * ey[0] + p[1] * ey[1] + p[2] * ey[2]];
    } else proj = p => [p[0], p[1]];
    const pos = o.coords || base;
    const all = [];
    for (const n of info.all_nodes) { all.push(proj(base[n])); if (o.coords && o.coords[n]) all.push(proj(o.coords[n])); }
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of all) { x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y); }
    const padL = d3 ? 16 : 44, padB = d3 ? 16 : 34, padR = 16, padT = top + 10;
    const sc = Math.min((W - padL - padR) / Math.max(x1 - x0, 1e-9), (H - padT - padB) / Math.max(y1 - y0, 1e-9));
    const cx = padL + (W - padL - padR - sc * (x1 - x0)) / 2, cy = padT + (H - padT - padB - sc * (y1 - y0)) / 2;
    const S = p => { const q = proj(p); return [cx + (q[0] - x0) * sc, cy + (y1 - q[1]) * sc]; };

    if (!d3) {          // light axes with ticks in metres
      const g = svgEl('g', {}, svg);
      const xt = ticks(x0, x1, Math.max(2, Math.round((x1 - x0) * sc / 80)));
      const yt = ticks(y0, y1, Math.max(2, Math.round((y1 - y0) * sc / 60)));
      const by = cy + (y1 - y0) * sc + 10;
      for (const v of xt.values) {
        if (v < x0 - 1e-9 || v > x1 + 1e-9) continue;
        const X = cx + (v - x0) * sc;
        svgEl('line', { x1: X, x2: X, y1: by, y2: by + 4, class: 'axis' }, g);
        text(g, X, by + 16, fmtTick(v, xt.step), { class: 'tick', 'text-anchor': 'middle' });
      }
      for (const v of yt.values) {
        if (v < y0 - 1e-9 || v > y1 + 1e-9) continue;
        const Yv = cy + (y1 - v) * sc;
        svgEl('line', { x1: cx - 14, x2: cx - 10, y1: Yv, y2: Yv, class: 'axis' }, g);
        text(g, cx - 17, Yv + 4, fmtTick(v, yt.step), { class: 'tick', 'text-anchor': 'end' });
      }
    } else {            // axis triad
      const g = svgEl('g', { class: 'triad' }, svg);
      const o0 = [36, H - 30];
      for (const [v, lab] of [[[1, 0, 0], 'X'], [[0, 1, 0], 'Y'], [[0, 0, 1], 'Z']]) {
        const q = proj(v), x = o0[0] + q[0] * 22, y = o0[1] - q[1] * 22;
        svgEl('line', { x1: o0[0], y1: o0[1], x2: x, y2: y, class: 'axis' }, g);
        text(g, x + (q[0] >= 0 ? 3 : -9), y + (q[1] >= 0 ? -2 : 10), lab, { class: 'tick' });
      }
    }

    const lw = d3 ? 1.6 : 2.4;
    const members = (P, cls, style) => {
      let d = '';
      for (const c of info.frames || info.columns.concat(info.beams)) {
        const a = S(P[c[1]]), b = S(P[c[2]]);
        d += `M${a[0].toFixed(1)} ${a[1].toFixed(1)}L${b[0].toFixed(1)} ${b[1].toFixed(1)}`;
      }
      svgEl('path', { d, class: cls, style }, svg);
    };
    const shells = (P, cls) => {
      for (const s of info.shells || []) svgEl('polygon', { points: s[1].map(n => S(P[n]).map(v => v.toFixed(1)).join(',')).join(' '), class: cls }, svg);
    };
    if (o.coords) {
      shells(base, 'ghost');
      shells(o.coords, 'shell');
      members(base, 'ghost', '');
      members(o.coords, 'member', `stroke:${o.color || 'var(--c1)'};stroke-width:${lw}`);
    } else {
      shells(base, 'shell');
      members(base, 'member', `stroke:${o.color || 'var(--c1)'};stroke-width:${lw}`);
    }

    // hinges: released member ends (open) and joint springs (filled), drawn a little inside the member
    for (const [n, other, kind] of info.hinges || []) {
      const pa = pos[n];
      let p = pa;
      if (other !== null && pos[other]) {
        const pb = pos[other], L = Math.hypot(...pb.map((v, i) => v - pa[i])) || 1, off = Math.min(0.04 * L, 0.35);
        p = pa.map((v, i) => v + (pb[i] - v) / L * off);
      }
      const s = S(p);
      svgEl('circle', { cx: s[0], cy: s[1], r: d3 ? 2.8 : 3.6, class: kind === 'pinned' ? 'hinge pinned' : 'hinge semi' }, svg);
    }
    for (const [n, sup] of info.supports || info.base_nodes.map(x => [x, null])) {
      const s = S(base[n]), r = d3 ? 3 : 5;
      const rot = sup && sup.slice(info.dim).some(Boolean), trans = !sup || sup.slice(0, info.dim).every(Boolean);
      if (!sup || (rot && trans)) svgEl('rect', { x: s[0] - r, y: s[1] - r, width: 2 * r, height: 2 * r, class: 'support' }, svg);
      else svgEl('polygon', { points: `${s[0]},${s[1]} ${s[0] - r - 1},${s[1] + 1.8 * r} ${s[0] + r + 1},${s[1] + 1.8 * r}`,
        class: trans ? 'support' : 'support open' }, svg);
    }
    if (d3) text(svg, W - 8, H - 8, 'drag to rotate · double-click to reset', { class: 'tick hint', 'text-anchor': 'end' });
  }

  // ------------------------------------------------------------------------
  // Cross-section drawing (mm): width horizontal, depth h vertical
  // ------------------------------------------------------------------------
  function section(el, s, title, lines) {
    el.innerHTML = '';
    el.classList.add('section-view');
    const k = s.kind;
    let bw, bh;
    const svg = svgEl('svg', { class: 'plot section-svg', role: 'img', 'aria-label': `${title} cross-section` });
    const g = svgEl('g', {}, svg);
    if (k === 'Steel I-section (fiber)') {
      const { d, bf, tw, tf } = s, w = bf / 2, hh = d / 2;
      const pts = [[-w, hh], [w, hh], [w, hh - tf], [tw / 2, hh - tf], [tw / 2, -hh + tf], [w, -hh + tf],
        [w, -hh], [-w, -hh], [-w, -hh + tf], [-tw / 2, -hh + tf], [-tw / 2, hh - tf], [-w, hh - tf]];
      svgEl('polygon', { points: pts.map(([x, y]) => `${x},${-y}`).join(' '), class: 'sec-steel' }, g);
      bw = bf; bh = d;
    } else {
      bw = s.b; bh = s.h;
      svgEl('rect', { x: -bw / 2, y: -bh / 2, width: bw, height: bh, class: 'sec-conc' }, g);
      if (k === 'RC rectangle (fiber)') {
        const c = s.cover;
        svgEl('rect', { x: -bw / 2 + c, y: -bh / 2 + c, width: bw - 2 * c, height: bh - 2 * c, class: 'sec-core' }, g);
        for (const [y, z] of FA.sections.rcBarLayout(s))
          svgEl('circle', { cx: z, cy: -y, r: s.d_bar / 2, class: 'sec-bar' }, g);
      }
    }
    const pad = 0.08 * Math.max(bw, bh), M = 0.6 * Math.max(bw, bh) + 2 * pad;
    svg.setAttribute('viewBox', `${-M} ${-M} ${2 * M} ${2 * M}`);
    const fs = 0.055 * 2 * M;
    const arrow = (x1, y1, x2, y2) => {
      svgEl('line', { x1, y1, x2, y2, class: 'dim', style: `stroke-width:${0.004 * M}` }, g);
      const a = Math.atan2(y2 - y1, x2 - x1), L = 0.035 * M;
      for (const [x, y, s2] of [[x1, y1, 1], [x2, y2, -1]]) {
        const p1 = [x + s2 * L * Math.cos(a - 0.4), y + s2 * L * Math.sin(a - 0.4)];
        const p2 = [x + s2 * L * Math.cos(a + 0.4), y + s2 * L * Math.sin(a + 0.4)];
        svgEl('polyline', { points: `${p1} ${x},${y} ${p2}`, class: 'dim', style: `stroke-width:${0.004 * M};fill:none` }, g);
      }
    };
    arrow(-bw / 2, bh / 2 + pad, bw / 2, bh / 2 + pad);
    text(g, 0, bh / 2 + 1.6 * pad + fs * 0.8, `${+bw.toFixed(1)}`, { class: 'dim-text', 'text-anchor': 'middle', style: `font-size:${fs}px` });
    arrow(-bw / 2 - pad, -bh / 2, -bw / 2 - pad, bh / 2);
    const t = text(g, 0, 0, `${+bh.toFixed(1)}`, { class: 'dim-text', 'text-anchor': 'middle', style: `font-size:${fs}px` });
    t.setAttribute('transform', `translate(${-bw / 2 - 1.6 * pad} 0) rotate(-90)`);
    const head = document.createElement('div');
    head.className = 'section-caption';
    head.innerHTML = `<strong>${title}</strong>` + lines.map(l => `<div>${l}</div>`).join('');
    el.appendChild(head);
    el.appendChild(svg);
  }

  FA.plots = { lineChart, frame, section, ticks, svgEl };
})(globalThis.FA = globalThis.FA || {});
