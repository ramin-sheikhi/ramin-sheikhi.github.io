/* Ground motion records: file parsing, sine wave and synthetic (filtered
 * white noise) accelerograms. Accelerations are returned in g.
 */
(function (FA) {
  'use strict';

  /** [acc array, dt]. fmt: 'PEER AT2', 'Single column', 'Time, Acc'. */
  function readRecord(text, fmt, dt) {
    const lines = String(text).trim().split(/\r?\n/);
    let vals, dtOut;
    if (fmt === 'PEER AT2') {
      if (lines.length < 5) throw new Error('could not read the PEER header (file has fewer than 5 lines)');
      const header = lines.slice(0, 4).join(' ');
      const mt = /NPTS\s*=\s*(\d+)\s*,?\s*DT\s*=\s*([\d.Ee+-]+)/i.exec(header);
      let npts;
      if (mt) { npts = parseInt(mt[1], 10); dtOut = parseFloat(mt[2]); }
      else {                         // older format: "5590   0.0050   NPTS, DT"
        const nums = lines[3].match(/[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?/g) || [];
        if (nums.length < 2) throw new Error('could not read NPTS and DT from the PEER header');
        npts = Math.trunc(parseFloat(nums[0])); dtOut = parseFloat(nums[1]);
      }
      vals = [];
      for (const ln of lines.slice(4)) for (const v of ln.trim().split(/\s+/)) if (v) vals.push(Number(v));
      vals = vals.slice(0, npts);
    } else {
      const rows = lines
        .filter(ln => ln.trim() && !ln.trimStart().startsWith('#') && /^\s*[-+.\d]/.test(ln))
        .map(ln => ln.trim().split(/[,\s;]+/).filter(Boolean).map(Number));
      if (fmt === 'Time, Acc') {
        const t = rows.map(r => r[0]);
        vals = rows.map(r => r[1]);
        if (t.length < 2) throw new Error('need at least two (time, acc) rows');
        dtOut = (t[t.length - 1] - t[0]) / (t.length - 1);
      } else {
        vals = rows.flat();
        dtOut = Number(dt);
      }
    }
    if (!vals.length || vals.some(v => !isFinite(v))) throw new Error('could not read numeric acceleration values');
    if (!(dtOut > 0)) throw new Error('could not read a positive time step dt');
    return [Float64Array.from(vals), dtOut];
  }

  function sineRecord(ampG, period, duration, dt) {
    const n = Math.ceil(duration / dt - 1e-9);
    const a = new Float64Array(n);
    for (let i = 0; i < n; i++) a[i] = ampG * Math.sin(2 * Math.PI * i * dt / period);
    return [a, dt];
  }

  /** In-place radix-2 FFT (re, im of length 2^k); inverse when inv is true. */
  function fft(re, im, inv) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) {
        let t = re[i]; re[i] = re[j]; re[j] = t;
        t = im[i]; im[i] = im[j]; im[j] = t;
      }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = 2 * Math.PI / len * (inv ? 1 : -1);
      const wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < len / 2; k++) {
          const a = i + k, b = a + len / 2;
          const xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr;
          re[b] = re[a] - xr; im[b] = im[a] - xi;
          re[a] += xr; im[a] += xi;
          const t = cr * wr - ci * wi;
          ci = cr * wi + ci * wr; cr = t;
        }
      }
    }
    if (inv) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  }

  /** Seeded normal random numbers (mulberry32 + Box-Muller). */
  function normalRng(seed) {
    let s = (seed >>> 0) || 1;
    const u = () => {
      s = (s + 0x6D2B79F5) >>> 0;
      let t = s;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    return () => {
      const a = u() || 1e-12, b = u();
      return Math.sqrt(-2 * Math.log(a)) * Math.cos(2 * Math.PI * b);
    };
  }

  /** Kanai-Tajimi filtered white noise with a trapezoidal envelope, scaled to PGA. */
  function syntheticRecord(pgaG, duration, dt, fg = 2.5, zg = 0.6, seed = 1) {
    const n = Math.trunc(duration / dt);
    if (n < 2) throw new Error('duration must be longer than two time steps');
    let N = 1;
    while (N < n) N <<= 1;
    const rnd = normalRng(seed);
    const re = new Float64Array(N), im = new Float64Array(N);
    for (let i = 0; i < N; i++) re[i] = rnd();
    fft(re, im, false);
    const wg = 2 * Math.PI * fg;
    for (let k = 0; k < N; k++) {
      const f = (k <= N / 2 ? k : N - k) / (N * dt);
      const w = 2 * Math.PI * f;
      const H2 = (wg ** 4 + 4 * zg ** 2 * wg ** 2 * w ** 2) /
                 ((wg ** 2 - w ** 2) ** 2 + 4 * zg ** 2 * wg ** 2 * w ** 2);
      const h = Math.sqrt(H2);
      re[k] *= h; im[k] *= h;
    }
    fft(re, im, true);
    const a = new Float64Array(n);
    const t1 = 0.1 * duration, t2 = 0.5 * duration;
    let peak = 0;
    for (let i = 0; i < n; i++) {
      const t = i * dt;
      const env = t < t1 ? (t / t1) ** 2 : t < t2 ? 1 : Math.exp(-3 * (t - t2) / (duration - t2));
      a[i] = re[i] * env;
      peak = Math.max(peak, Math.abs(a[i]));
    }
    for (let i = 0; i < n; i++) a[i] *= pgaG / peak;
    return [a, dt];
  }

  FA.records = { readRecord, sineRecord, syntheticRecord, fft };
})(globalThis.FA = globalThis.FA || {});
