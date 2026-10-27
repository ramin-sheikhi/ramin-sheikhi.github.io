/* Frame elements in 2D (3 dof/node) and 3D (6 dof/node).
 *
 *   beamIntegration(type, n)   Lobatto, Legendre, Radau, NewtonCotes on [0, 1]
 *   Transf                     Linear / PDelta coordinate transformation
 *   ElasticBeam                elasticBeamColumn
 *   DispBeam                   dispBeamColumn (cubic Hermite, sections at IPs)
 *   ForceBeam                  forceBeamColumn (flexibility, iterative state
 *                              determination with subdivision)
 *   ZeroLengthRot              rotational spring between two coincident nodes
 *
 * Basic system (OpenSees): 2D q = [N, Mi, Mj]; 3D q = [N, Mzi, Mzj, Myi, Myj, T].
 * Every element exposes update(ug, lam) -> bool, K(), F(), Kinit(), commit(),
 * revert(). ug is the element displacement vector in global axes, lam the
 * factor applied to its distributed (gravity) load.
 */
(function (FA) {
  'use strict';
  const { inv } = FA.linalg;

  // ------------------------------------------------------------------------
  // Integration rules on [0, 1]
  // ------------------------------------------------------------------------
  const ruleCache = {};

  function legendre(n, x) {          // [P_n(x), P_n'(x)]
    if (n === 0) return [1, 0];
    let p0 = 1, p1 = x;
    for (let k = 2; k <= n; k++) {
      const p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
      p0 = p1; p1 = p2;
    }
    const dp = Math.abs(1 - x * x) < 1e-14 ? 0.5 * n * (n + 1) * (x > 0 ? 1 : (n % 2 ? 1 : -1)) : n * (p0 - x * p1) / (1 - x * x);
    return [p1, dp];
  }

  /** Roots of f in the open interval (-1, 1) by scanning + bisection. */
  function roots(f) {
    const out = [], N = 4000;
    let xa = -1 + 1e-9, fa = f(xa);
    for (let i = 1; i <= N; i++) {
      const xb = i === N ? 1 - 1e-9 : -1 + 2 * i / N, fb = f(xb);
      if (fa === 0) out.push(xa);
      else if (fa * fb < 0) {
        let a = xa, b = xb, va = fa;
        for (let k = 0; k < 80; k++) {
          const m = 0.5 * (a + b), vm = f(m);
          if (va * vm <= 0) b = m; else { a = m; va = vm; }
        }
        out.push(0.5 * (a + b));
      }
      xa = xb; fa = fb;
    }
    return out;
  }

  function beamIntegration(type, n) {
    const key = type + n;
    if (ruleCache[key]) return ruleCache[key];
    let x;                                   // points on [-1, 1]
    if (type === 'Legendre') x = roots(t => legendre(n, t)[0]);
    else if (type === 'Lobatto') x = [-1].concat(n > 2 ? roots(t => legendre(n - 1, t)[1]) : [], [1]);
    else if (type === 'Radau') x = [-1].concat(roots(t => legendre(n - 1, t)[0] + legendre(n, t)[0]));
    else if (type === 'NewtonCotes') { x = []; for (let i = 0; i < n; i++) x.push(-1 + 2 * i / (n - 1)); }
    else throw new Error('Unknown beam integration ' + type);
    if (x.length !== n) throw new Error(`${type} integration: found ${x.length} of ${n} points`);
    const xi = x.map(t => 0.5 * (t + 1));
    // interpolatory weights on [0, 1]: Σ w_i xi_i^k = 1/(k+1)
    const V = new Float64Array(n * n);
    for (let k = 0; k < n; k++) for (let i = 0; i < n; i++) V[k * n + i] = Math.pow(xi[i], k);
    const Vi = inv(V, n);
    const w = [];
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += Vi[i * n + k] / (k + 1);
      w.push(s);
    }
    return (ruleCache[key] = { xi, w });
  }

  // ------------------------------------------------------------------------
  // Coordinate transformation
  // ------------------------------------------------------------------------
  function cross(a, b) {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  }
  function unit(a) {
    const l = Math.hypot(...a);
    return a.map(v => v / l);
  }

  class Transf {
    constructor(dim, xi, xj, vecxz, pdelta) {
      this.dim = dim;
      const ndf = dim === 2 ? 3 : 6, n = 2 * ndf, nb = dim === 2 ? 3 : 6;
      Object.assign(this, { ndf, n, nb, pdelta: !!pdelta });
      const d = xj.map((v, i) => v - xi[i]);
      const L = Math.hypot(...d);
      if (!(L > 0)) throw new Error('Element of zero length');
      this.L = L;
      const T = new Float64Array(n * n);
      if (dim === 2) {
        const c = d[0] / L, s = d[1] / L;
        for (const o of [0, 3]) {
          T[o * n + o] = c; T[o * n + o + 1] = s;
          T[(o + 1) * n + o] = -s; T[(o + 1) * n + o + 1] = c;
          T[(o + 2) * n + o + 2] = 1;
        }
      } else {
        const ex = d.map(v => v / L);
        const ey = unit(cross(vecxz, ex));
        const ez = cross(ex, ey);
        const R = [ex, ey, ez];
        for (const o of [0, 3, 6, 9])
          for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) T[(o + r) * n + o + c] = R[r][c];
        this.R = R;
      }
      this.T = T;
      // A: basic deformations from local displacements
      const A = new Float64Array(nb * n), iL = 1 / L;
      const set = (r, c, v) => { A[r * n + c] = v; };
      if (dim === 2) {
        set(0, 0, -1); set(0, 3, 1);
        set(1, 1, iL); set(1, 4, -iL); set(1, 2, 1);
        set(2, 1, iL); set(2, 4, -iL); set(2, 5, 1);
      } else {
        set(0, 0, -1); set(0, 6, 1);
        set(1, 1, iL); set(1, 7, -iL); set(1, 5, 1);
        set(2, 1, iL); set(2, 7, -iL); set(2, 11, 1);
        set(3, 2, -iL); set(3, 8, iL); set(3, 4, 1);
        set(4, 2, -iL); set(4, 8, iL); set(4, 10, 1);
        set(5, 3, -1); set(5, 9, 1);
      }
      const AT = new Float64Array(nb * n);
      for (let r = 0; r < nb; r++)
        for (let c = 0; c < n; c++) {
          let s = 0;
          for (let k = 0; k < n; k++) s += A[r * n + k] * T[k * n + c];
          AT[r * n + c] = s;
        }
      this.AT = AT;
      // P-Delta: difference of the transverse local rows (y, and z in 3D)
      const g = (pi, pj) => {
        const v = new Float64Array(n);
        for (let c = 0; c < n; c++) v[c] = T[pi * n + c] - T[pj * n + c];
        return v;
      };
      this.g = dim === 2 ? [g(1, 4)] : [g(1, 7), g(2, 8)];
    }

    basic(ug, v) {
      const { nb, n, AT } = this;
      v = v || new Float64Array(nb);
      for (let r = 0; r < nb; r++) {
        let s = 0;
        for (let c = 0; c < n; c++) s += AT[r * n + c] * ug[c];
        v[r] = s;
      }
      return v;
    }

    /** Global forces from basic forces q and the (global) end-load vector p0g * lam. */
    force(q, ug, p0g, lam) {
      const { nb, n, AT } = this;
      const F = new Float64Array(n);
      for (let r = 0; r < nb; r++) {
        const qr = q[r];
        if (qr === 0) continue;
        for (let c = 0; c < n; c++) F[c] += AT[r * n + c] * qr;
      }
      if (p0g && lam) for (let c = 0; c < n; c++) F[c] += lam * p0g[c];
      if (this.pdelta) {
        const NL = q[0] / this.L;
        for (const g of this.g) {
          let gu = 0;
          for (let c = 0; c < n; c++) gu += g[c] * ug[c];
          for (let c = 0; c < n; c++) F[c] += NL * gu * g[c];
        }
      }
      return F;
    }

    stiff(kb, q) {
      const { nb, n, AT } = this;
      const tmp = new Float64Array(nb * n);   // kb * AT
      for (let r = 0; r < nb; r++)
        for (let k = 0; k < nb; k++) {
          const v = kb[r * nb + k];
          if (v === 0) continue;
          for (let c = 0; c < n; c++) tmp[r * n + c] += v * AT[k * n + c];
        }
      const K = new Float64Array(n * n);
      for (let r = 0; r < nb; r++)
        for (let i = 0; i < n; i++) {
          const a = AT[r * n + i];
          if (a === 0) continue;
          for (let j = 0; j < n; j++) K[i * n + j] += a * tmp[r * n + j];
        }
      if (this.pdelta && q) {
        const NL = q[0] / this.L;
        for (const g of this.g)
          for (let i = 0; i < n; i++) {
            if (g[i] === 0) continue;
            for (let j = 0; j < n; j++) K[i * n + j] += NL * g[i] * g[j];
          }
      }
      return K;
    }

    /** Global end-load vector of a uniform load wy (local y). */
    uniformP0(wy) {
      const { n, T, L } = this;
      const pl = new Float64Array(n);
      const V = -0.5 * wy * L;
      if (this.dim === 2) { pl[1] = V; pl[4] = V; } else { pl[1] = V; pl[7] = V; }
      const p = new Float64Array(n);
      for (let c = 0; c < n; c++) {
        let s = 0;
        for (let r = 0; r < n; r++) s += T[r * n + c] * pl[r];
        p[c] = s;
      }
      return p;
    }
  }

  // ------------------------------------------------------------------------
  // Beam elements
  // ------------------------------------------------------------------------
  class BeamBase {
    constructor(transf, wy) {
      this.tr = transf;
      this.n = transf.n;
      this.nb = transf.nb;
      this.wy = wy || 0;
      this.p0g = this.wy ? transf.uniformP0(this.wy) : null;
      this.ug = new Float64Array(this.n);
      this.Cug = new Float64Array(this.n);
      this.lam = 0;
      this.Clam = 0;
      this.q = new Float64Array(this.nb);
    }
    /** fixed-end basic forces of the uniform load for lam = 1 */
    q0ref() {
      const q0 = new Float64Array(this.nb);
      if (this.wy) {
        const M = this.wy * this.tr.L ** 2 / 12;
        q0[1] = -M; q0[2] = M;
      }
      return q0;
    }
    F() { return this.tr.force(this.q, this.ug, this.p0g, this.lam); }
    K() { return this.tr.stiff(this.kb, this.q); }
    Kinit() { return this.tr.stiff(this.kb0, null); }
  }

  function elasticKb(dim, sp, L) {
    const nb = dim === 2 ? 3 : 6;
    const kb = new Float64Array(nb * nb);
    kb[0] = sp.E * sp.A / L;
    const bend = (i, EI) => {
      kb[i * nb + i] = 4 * EI / L; kb[(i + 1) * nb + i + 1] = 4 * EI / L;
      kb[i * nb + i + 1] = 2 * EI / L; kb[(i + 1) * nb + i] = 2 * EI / L;
    };
    bend(1, sp.E * sp.Iz);
    if (dim === 3) { bend(3, sp.E * sp.Iy); kb[35] = sp.G * sp.J / L; }
    return kb;
  }

  class ElasticBeam extends BeamBase {
    constructor(dim, transf, sp, wy) {
      super(transf, wy);
      this.kb = elasticKb(dim, sp, transf.L);
      this.kb0 = this.kb;
      this.q0 = this.q0ref();
      this.v = new Float64Array(this.nb);
    }
    update(ug, lam) {
      this.ug.set(ug);
      this.lam = lam;
      const { nb, kb, q, q0 } = this;
      const v = this.tr.basic(ug, this.v);
      for (let r = 0; r < nb; r++) {
        let s = lam * q0[r];
        for (let c = 0; c < nb; c++) s += kb[r * nb + c] * v[c];
        q[r] = s;
      }
      return true;
    }
    commit() { this.Cug.set(this.ug); this.Clam = this.lam; }
    revert() { this.update(this.Cug, this.Clam); }
  }

  /** Section deformation-interpolation matrix rows of a displacement-based element. */
  function dispB(dim, xi, L) {
    const iL = 1 / L, a = (6 * xi - 4) * iL, b = (6 * xi - 2) * iL;
    if (dim === 2) return Float64Array.from([iL, 0, 0, 0, a, b]);                // 2 x 3
    const B = new Float64Array(24);                                              // 4 x 6
    B[0] = iL; B[6 + 1] = a; B[6 + 2] = b; B[12 + 3] = a; B[12 + 4] = b; B[18 + 5] = iL;
    return B;
  }

  /** Force-interpolation matrix rows b(ξ) of a force-based element. */
  function forceB(dim, xi) {
    if (dim === 2) return Float64Array.from([1, 0, 0, 0, xi - 1, xi]);
    const b = new Float64Array(24);
    b[0] = 1; b[6 + 1] = xi - 1; b[6 + 2] = xi; b[12 + 3] = xi - 1; b[12 + 4] = xi; b[18 + 5] = 1;
    return b;
  }

  class DispBeam extends BeamBase {
    constructor(dim, transf, makeSection, rule, wy) {
      super(transf, wy);
      const L = transf.L;
      this.secs = rule.xi.map(() => makeSection());
      this.ord = this.secs[0].order;
      this.B = rule.xi.map(x => dispB(dim, x, L));
      this.wL = rule.w.map(w => w * L);
      this.q0 = this.q0ref();
      this.v = new Float64Array(this.nb);
      this.e = new Float64Array(this.ord);
      this.kb = new Float64Array(this.nb * this.nb);
      this.integrate(0);
      this.kb0 = Float64Array.from(this.kb);
    }
    integrate(lam) {
      const { nb, ord, q, kb, q0 } = this;
      for (let r = 0; r < nb; r++) q[r] = lam * q0[r];
      kb.fill(0);
      this.secs.forEach((sec, i) => {
        const B = this.B[i], wL = this.wL[i], s = sec.s, ks = sec.ks;
        for (let a = 0; a < ord; a++) {
          const sa = s[a] * wL;
          for (let r = 0; r < nb; r++) q[r] += B[a * nb + r] * sa;
        }
        for (let a = 0; a < ord; a++)
          for (let b = 0; b < ord; b++) {
            const k = ks[a * ord + b] * wL;
            if (k === 0) continue;
            for (let r = 0; r < nb; r++) {
              const Bar = B[a * nb + r];
              if (Bar === 0) continue;
              for (let c = 0; c < nb; c++) kb[r * nb + c] += Bar * k * B[b * nb + c];
            }
          }
      });
    }
    update(ug, lam) {
      this.ug.set(ug);
      this.lam = lam;
      const { nb, ord, e } = this;
      const v = this.tr.basic(ug, this.v);
      this.secs.forEach((sec, i) => {
        const B = this.B[i];
        for (let a = 0; a < ord; a++) {
          let s = 0;
          for (let r = 0; r < nb; r++) s += B[a * nb + r] * v[r];
          e[a] = s;
        }
        sec.setTrial(e);
      });
      this.integrate(lam);
      return true;
    }
    commit() {
      this.Cug.set(this.ug); this.Clam = this.lam;
      for (const s of this.secs) s.commit();
    }
    revert() {
      this.ug.set(this.Cug); this.lam = this.Clam;
      for (const s of this.secs) s.revert();
      this.integrate(this.lam);
    }
  }

  class ForceBeam extends BeamBase {
    constructor(dim, transf, makeSection, rule, wy) {
      super(transf, wy);
      const L = transf.L;
      this.dim = dim;
      this.secs = rule.xi.map(() => makeSection());
      this.ord = this.secs[0].order;
      this.b = rule.xi.map(x => forceB(dim, x));
      this.wL = rule.w.map(w => w * L);
      // section forces of the uniform load (simply supported basic system), lam = 1
      this.sp = rule.xi.map(x => {
        const s = new Float64Array(this.ord);
        const xx = x * L;
        if (this.wy) s[1] = this.wy * 0.5 * xx * (xx - L);
        return s;
      });
      const ns = this.secs.length, ord = this.ord, nb = this.nb;
      this.es = this.secs.map(() => new Float64Array(ord));
      this.sr = this.secs.map(() => new Float64Array(ord));
      this.fs = this.secs.map(sec => inv(sec.ks, ord));
      if (this.fs.some(f => !f)) throw new Error('forceBeamColumn: singular initial section stiffness');
      this.fs0 = this.fs.map(f => Float64Array.from(f));
      this.kv = this.flexToStiff();
      if (!this.kv) throw new Error('forceBeamColumn: singular initial element flexibility');
      this.kb0 = Float64Array.from(this.kv);
      this.vCur = new Float64Array(nb);
      this.v = new Float64Array(nb);
      this.C = this.snapshot();
      this.tol = 1e-12;
      this.maxIter = 10;
      this.ns = ns;
    }
    get kb() { return this.kv; }
    flexToStiff() {
      const nb = this.nb, ord = this.ord;
      const f = new Float64Array(nb * nb);
      this.secs.forEach((sec, i) => this.addFlex(f, this.b[i], this.fs[i], this.wL[i]));
      return inv(f, nb);
    }
    addFlex(f, b, fs, wL) {
      const nb = this.nb, ord = this.ord;
      for (let a = 0; a < ord; a++)
        for (let c = 0; c < ord; c++) {
          const v = fs[a * ord + c] * wL;
          if (v === 0) continue;
          for (let r = 0; r < nb; r++) {
            const bar = b[a * nb + r];
            if (bar === 0) continue;
            for (let k = 0; k < nb; k++) f[r * nb + k] += bar * v * b[c * nb + k];
          }
        }
    }
    snapshot() {
      return {
        q: Float64Array.from(this.q), kv: Float64Array.from(this.kv), vCur: Float64Array.from(this.vCur || new Float64Array(this.nb)),
        es: this.es.map(a => Float64Array.from(a)), sr: this.sr.map(a => Float64Array.from(a)),
        fs: this.fs.map(a => Float64Array.from(a)), lam: this.lam,
      };
    }
    restore(S) {
      this.q.set(S.q); this.kv.set(S.kv); this.vCur.set(S.vCur); this.lam = S.lam;
      S.es.forEach((a, i) => this.es[i].set(a));
      S.sr.forEach((a, i) => this.sr[i].set(a));
      S.fs.forEach((a, i) => this.fs[i].set(a));
    }
    update(ug, lam) {
      this.ug.set(ug);
      const nb = this.nb;
      const v = this.tr.basic(ug, this.v);
      let same = lam === this.lam;
      for (let r = 0; r < nb && same; r++) if (v[r] !== this.vCur[r]) same = false;
      if (same) return true;
      const v0 = Float64Array.from(this.vCur), lam0 = this.lam;
      const dv = new Float64Array(nb);
      for (let r = 0; r < nb; r++) dv[r] = v[r] - v0[r];
      let done = 0, step = 1;
      const vt = new Float64Array(nb);
      while (done < 1 - 1e-12) {
        const to = Math.min(1, done + step);
        for (let r = 0; r < nb; r++) vt[r] = v0[r] + to * dv[r];
        const lt = lam0 + to * (lam - lam0);
        const S = this.snapshot();
        // as OpenSees: Newton, then initial-flexibility iterations, then
        // initial flexibility on the first iteration only; then subdivide
        let ok = false;
        for (let l = 0; l < 3 && !ok; l++) {
          if (l) this.restore(S);
          ok = this.iterate(vt, lt, l);
        }
        if (ok) {
          done = to;
        } else {
          this.restore(S);
          step /= 10;
          if (step < 1e-4) {
            this.failMsg = `WARNING forceBeamColumn: failed to get compatible element forces & deformations`;
            return false;
          }
        }
      }
      return true;
    }
    /** l = 0: Newton; 1: initial flexibility; 2: initial flexibility on the first iteration. */
    iterate(vt, lt, l) {
      const { nb, ord, secs } = this;
      const dq = new Float64Array(nb), dvv = new Float64Array(nb);
      for (let r = 0; r < nb; r++) dvv[r] = vt[r] - this.vCur[r];
      mv(l ? this.kb0 : this.kv, dvv, dq, nb);
      const s = new Float64Array(ord), ds = new Float64Array(ord), res = new Float64Array(ord);
      const f = new Float64Array(nb * nb), vr = new Float64Array(nb);
      const numIter = l === 1 ? 10 * this.maxIter : this.maxIter;
      for (let it = 0; it < numIter; it++) {
        const useInit = l === 1 || (l === 2 && it === 0);
        for (let r = 0; r < nb; r++) this.q[r] += dq[r];
        f.fill(0); vr.fill(0);
        for (let i = 0; i < secs.length; i++) {
          const b = this.b[i], es = this.es[i], sr = this.sr[i], sp = this.sp[i], wL = this.wL[i];
          const fsCur = this.fs[i], fs = useInit ? this.fs0[i] : fsCur;
          for (let a = 0; a < ord; a++) {
            let t = lt * sp[a];
            for (let r = 0; r < nb; r++) t += b[a * nb + r] * this.q[r];
            s[a] = t;
            ds[a] = t - sr[a];
          }
          for (let a = 0; a < ord; a++) {
            let t = 0;
            for (let c = 0; c < ord; c++) t += fs[a * ord + c] * ds[c];
            es[a] += t;
          }
          secs[i].setTrial(es);
          sr.set(secs[i].s);
          if (!inv(secs[i].ks, ord, fsCur)) return false;
          for (let a = 0; a < ord; a++) ds[a] = s[a] - sr[a];
          for (let a = 0; a < ord; a++) {
            let t = 0;
            for (let c = 0; c < ord; c++) t += fs[a * ord + c] * ds[c];
            res[a] = es[a] + t;
          }
          this.addFlex(f, b, fsCur, wL);
          for (let a = 0; a < ord; a++) {
            const t = res[a] * wL;
            if (t === 0) continue;
            for (let r = 0; r < nb; r++) vr[r] += b[a * nb + r] * t;
          }
        }
        if (!inv(f, nb, this.kv)) return false;
        for (let r = 0; r < nb; r++) dvv[r] = vt[r] - vr[r];
        mv(useInit ? this.kb0 : this.kv, dvv, dq, nb);
        let dW = 0;
        for (let r = 0; r < nb; r++) dW += dq[r] * dvv[r];
        if (!isFinite(dW)) return false;
        if (Math.abs(dW) < this.tol) {
          this.vCur.set(vt);
          this.lam = lt;
          return true;
        }
      }
      return false;
    }
    commit() {
      this.Cug.set(this.ug);
      for (const s of this.secs) s.commit();
      this.C = this.snapshot();
    }
    revert() {
      this.ug.set(this.Cug);
      for (const s of this.secs) s.revert();
      this.restore(this.C);
    }
  }

  function mv(A, x, y, n) {
    for (let r = 0; r < n; r++) {
      let s = 0;
      for (let c = 0; c < n; c++) s += A[r * n + c] * x[c];
      y[r] = s;
    }
  }

  /**
   * zeroLength element between nodes a and b: one uniaxial material per
   * spring, acting on the relative displacement of dof `dof` (0-based, global
   * axes). springs: [{dof, make}] where make() returns a new material.
   */
  class ZeroLength {
    constructor(ndf, springs) {
      this.ndf = ndf;
      this.n = 2 * ndf;
      this.springs = springs.map(s => ({ dof: s.dof, mat: s.make() }));
      this.noRayleigh = true;          // zeroLength elements carry no Rayleigh damping by default
      this.ug = new Float64Array(this.n);
      this.Cug = new Float64Array(this.n);
    }
    update(ug) {
      this.ug.set(ug);
      for (const s of this.springs) s.mat.setTrial(ug[this.ndf + s.dof] - ug[s.dof]);
      return true;
    }
    kmat(tangent) {
      const n = this.n, K = new Float64Array(n * n);
      for (const s of this.springs) {
        const k = tangent ? s.mat.tan : s.mat.E0, i = s.dof, j = this.ndf + s.dof;
        K[i * n + i] += k; K[j * n + j] += k; K[i * n + j] -= k; K[j * n + i] -= k;
      }
      return K;
    }
    K() { return this.kmat(true); }
    Kinit() { return this.kmat(false); }
    F() {
      const F = new Float64Array(this.n);
      for (const s of this.springs) {
        F[s.dof] -= s.mat.sig;
        F[this.ndf + s.dof] += s.mat.sig;
      }
      return F;
    }
    commit() { for (const s of this.springs) s.mat.commit(); this.Cug.set(this.ug); }
    revert() { for (const s of this.springs) s.mat.revert(); this.ug.set(this.Cug); }
  }

  /** Rotational spring on dof `dof` (0-based within a node) between nodes a and b. */
  class ZeroLengthRot extends ZeroLength {
    constructor(ndf, dof, makeMat) { super(ndf, [{ dof, make: makeMat }]); }
  }

  // ------------------------------------------------------------------------
  // Linear elastic area elements (constant stiffness in global axes)
  // ------------------------------------------------------------------------
  class LinearElement {
    constructor(Kg, n) {
      this.Kg = Kg;
      this.n = n;
      this.ug = new Float64Array(n);
      this.Cug = new Float64Array(n);
    }
    update(ug) { this.ug.set(ug); return true; }
    K() { return this.Kg; }
    Kinit() { return this.Kg; }
    F() {
      const n = this.n, F = new Float64Array(n), K = this.Kg, u = this.ug;
      for (let r = 0; r < n; r++) {
        let s = 0;
        for (let c = 0; c < n; c++) s += K[r * n + c] * u[c];
        F[r] = s;
      }
      return F;
    }
    commit() { this.Cug.set(this.ug); }
    revert() { this.ug.set(this.Cug); }
  }

  const QXI = [-1, 1, 1, -1], QETA = [-1, -1, 1, 1];
  const GP = 1 / Math.sqrt(3);
  const GAUSS4 = [[-GP, -GP], [GP, -GP], [GP, GP], [-GP, GP]];

  /** Bilinear shape functions, derivatives, Jacobian and its inverse at (ξ, η). */
  function quadShape(xy, xi, eta) {
    const N = [], Nxi = [], Neta = [];
    for (let i = 0; i < 4; i++) {
      N.push(0.25 * (1 + QXI[i] * xi) * (1 + QETA[i] * eta));
      Nxi.push(0.25 * QXI[i] * (1 + QETA[i] * eta));
      Neta.push(0.25 * QETA[i] * (1 + QXI[i] * xi));
    }
    let j00 = 0, j01 = 0, j10 = 0, j11 = 0;
    for (let i = 0; i < 4; i++) {
      j00 += Nxi[i] * xy[i][0]; j01 += Nxi[i] * xy[i][1];
      j10 += Neta[i] * xy[i][0]; j11 += Neta[i] * xy[i][1];
    }
    const det = j00 * j11 - j01 * j10;
    const inv = [j11 / det, -j01 / det, -j10 / det, j00 / det];
    const Nx = [], Ny = [];
    for (let i = 0; i < 4; i++) {
      Nx.push(inv[0] * Nxi[i] + inv[1] * Neta[i]);
      Ny.push(inv[2] * Nxi[i] + inv[3] * Neta[i]);
    }
    return { N, Nxi, Neta, J: [j00, j01, j10, j11], det, inv, Nx, Ny };
  }

  function quadArea(xy) {
    let a = 0;
    for (let i = 0; i < 4; i++) { const p = xy[i], q = xy[(i + 1) % 4]; a += p[0] * q[1] - q[0] * p[1]; }
    return 0.5 * a;
  }

  /** k += Bᵀ D B · w for B (m x n) and D (m x m), all row-major. */
  function addBtDB(k, B, D, m, n, w) {
    const DB = new Float64Array(m * n);
    for (let a = 0; a < m; a++)
      for (let b = 0; b < m; b++) {
        const d = D[a * m + b];
        if (d === 0) continue;
        for (let c = 0; c < n; c++) DB[a * n + c] += d * B[b * n + c];
      }
    for (let a = 0; a < m; a++)
      for (let r = 0; r < n; r++) {
        const br = B[a * n + r] * w;
        if (br === 0) continue;
        for (let c = 0; c < n; c++) k[r * n + c] += br * DB[a * n + c];
      }
  }

  /** Plane-stress constitutive matrix (3 x 3) times thickness t. */
  function planeStressD(E, nu, t) {
    const c = E * t / (1 - nu * nu);
    return Float64Array.from([c, c * nu, 0, c * nu, c, 0, 0, 0, c * (1 - nu) / 2]);
  }

  /**
   * Membrane stiffness (8 x 8, order u1 v1 ... u4 v4) of a 4-node quad with
   * incompatible modes (Wilson-Taylor Q6), statically condensed. xy are
   * counter-clockwise local coordinates.
   */
  function q6Membrane(xy, D) {
    const c0 = quadShape(xy, 0, 0);
    const Kcc = new Float64Array(64), Kci = new Float64Array(32), Kii = new Float64Array(16);
    for (const [xi, eta] of GAUSS4) {
      const s = quadShape(xy, xi, eta);
      const Bc = new Float64Array(24);
      for (let i = 0; i < 4; i++) {
        Bc[2 * i] = s.Nx[i]; Bc[8 + 2 * i + 1] = s.Ny[i];
        Bc[16 + 2 * i] = s.Ny[i]; Bc[16 + 2 * i + 1] = s.Nx[i];
      }
      // incompatible modes 1 - ξ², 1 - η² with the centre Jacobian (Taylor's correction)
      const r = c0.det / s.det, dP = [[-2 * xi, 0], [0, -2 * eta]];
      const Px = dP.map(([a, b]) => (c0.inv[0] * a + c0.inv[1] * b) * r);
      const Py = dP.map(([a, b]) => (c0.inv[2] * a + c0.inv[3] * b) * r);
      const Bi = Float64Array.from([Px[0], Px[1], 0, 0, 0, 0, Py[0], Py[1], Py[0], Py[1], Px[0], Px[1]]);
      const w = s.det;
      addBtDB(Kcc, Bc, D, 3, 8, w);
      addBtDB(Kii, Bi, D, 3, 4, w);
      const DBi = new Float64Array(12);
      for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) for (let c = 0; c < 4; c++) DBi[a * 4 + c] += D[a * 3 + b] * Bi[b * 4 + c];
      for (let a = 0; a < 3; a++) for (let r2 = 0; r2 < 8; r2++) for (let c = 0; c < 4; c++) Kci[r2 * 4 + c] += Bc[a * 8 + r2] * DBi[a * 4 + c] * w;
    }
    const Kinv = inv(Kii, 4);
    if (!Kinv) throw new Error('Shell element: singular incompatible-mode matrix (distorted element?)');
    for (let r = 0; r < 8; r++)
      for (let c = 0; c < 8; c++) {
        let s = 0;
        for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) s += Kci[r * 4 + a] * Kinv[a * 4 + b] * Kci[c * 4 + b];
        Kcc[r * 8 + c] -= s;
      }
    return Kcc;
  }

  /**
   * MITC4 plate bending stiffness (12 x 12, order w θx θy per node) with
   * the assumed transverse shear strain field of Bathe and Dvorkin.
   * Rotations are right-handed about the local axes: βx = θy, βy = −θx.
   */
  function mitc4Bending(xy, E, nu, t) {
    const Db = planeStressD(E, nu, 1).map(v => v * t ** 3 / 12);
    const Ds = 5 / 6 * E / (2 * (1 + nu)) * t;
    // covariant transverse shear rows γξz, γηz at (ξ, η)
    const cov = (xi, eta) => {
      const s = quadShape(xy, xi, eta), [xx, xy_, ex, ey] = s.J;
      const rx = new Float64Array(12), re = new Float64Array(12);
      for (let i = 0; i < 4; i++) {
        rx[3 * i] = s.Nxi[i]; rx[3 * i + 2] = s.N[i] * xx; rx[3 * i + 1] = -s.N[i] * xy_;
        re[3 * i] = s.Neta[i]; re[3 * i + 2] = s.N[i] * ex; re[3 * i + 1] = -s.N[i] * ey;
      }
      return [rx, re];
    };
    const A = cov(0, 1)[0], C = cov(0, -1)[0], B = cov(-1, 0)[1], Dd = cov(1, 0)[1];
    const K = new Float64Array(144);
    for (const [xi, eta] of GAUSS4) {
      const s = quadShape(xy, xi, eta);
      const Bb = new Float64Array(36);
      for (let i = 0; i < 4; i++) {
        Bb[3 * i + 2] = s.Nx[i];                          // κx = ∂βx/∂x
        Bb[12 + 3 * i + 1] = -s.Ny[i];                    // κy = ∂βy/∂y
        Bb[24 + 3 * i + 2] = s.Ny[i]; Bb[24 + 3 * i + 1] = -s.Nx[i];
      }
      addBtDB(K, Bb, Db, 3, 12, s.det);
      const Bs = new Float64Array(24);
      for (let c = 0; c < 12; c++) {
        const gx = 0.5 * (1 + eta) * A[c] + 0.5 * (1 - eta) * C[c];
        const ge = 0.5 * (1 + xi) * Dd[c] + 0.5 * (1 - xi) * B[c];
        Bs[c] = s.inv[0] * gx + s.inv[1] * ge;
        Bs[12 + c] = s.inv[2] * gx + s.inv[3] * ge;
      }
      addBtDB(K, Bs, Float64Array.from([Ds, 0, 0, Ds]), 2, 12, s.det);
    }
    return K;
  }

  /** Counter-clockwise reordering of 4 planar points; returns the permutation. */
  function ccw(xy) {
    return quadArea(xy) >= 0 ? [0, 1, 2, 3] : [0, 3, 2, 1];
  }

  /**
   * 2D plane-stress quad (Q6) for frames in the X-Y plane. Nodes carry
   * 3 dofs (ux, uy, rz); the rotation is not used by the element.
   */
  function quad2D(coords, E, nu, t) {
    const p = ccw(coords), xy = p.map(i => coords[i]);
    if (!(quadArea(xy) > 0)) throw new Error('Shell element has zero area');
    const km = q6Membrane(xy, planeStressD(E, nu, t));
    const n = 12, K = new Float64Array(n * n);
    for (let a = 0; a < 8; a++)
      for (let b = 0; b < 8; b++) {
        const r = 3 * p[a >> 1] + (a & 1), c = 3 * p[b >> 1] + (b & 1);
        K[r * n + c] += km[a * 8 + b];
      }
    return new LinearElement(K, n);
  }

  /**
   * Flat 4-node shell for 3D frames (6 dofs per node): Q6 membrane, MITC4
   * plate bending and a Hughes-Brezzi drilling stiffness. Mildly warped
   * elements are projected on their mean plane.
   */
  function shell3D(coords, E, nu, t, drill = 1e-3) {
    const c = [0, 1, 2].map(k => 0.25 * coords.reduce((s, p) => s + p[k], 0));
    const sub = (a, b) => a.map((v, i) => v - b[i]);
    let nrm = cross(sub(coords[2], coords[0]), sub(coords[3], coords[1]));
    if (!(Math.hypot(...nrm) > 0)) throw new Error('Shell element has zero area');
    nrm = unit(nrm);
    const mid = (i, j) => coords[i].map((v, k) => 0.5 * (v + coords[j][k]));
    let e1 = sub(mid(1, 2), mid(0, 3));
    const dn = e1.reduce((s, v, k) => s + v * nrm[k], 0);
    e1 = unit(e1.map((v, k) => v - dn * nrm[k]));
    const e2 = cross(nrm, e1);
    const R = [e1, e2, nrm];
    const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const xy = coords.map(p => { const d = sub(p, c); return [dot(d, e1), dot(d, e2)]; });
    const warp = Math.max(...coords.map(p => Math.abs(dot(sub(p, c), nrm))));
    const A = quadArea(xy);
    if (!(A > 0)) throw new Error('Shell element has zero area');
    const km = q6Membrane(xy, planeStressD(E, nu, t));
    const kb = mitc4Bending(xy, E, nu, t);
    const Kl = new Float64Array(576);
    const add = (r, cc, v) => { Kl[r * 24 + cc] += v; };
    for (let a = 0; a < 8; a++) for (let b = 0; b < 8; b++) add(6 * (a >> 1) + (a & 1), 6 * (b >> 1) + (b & 1), km[a * 8 + b]);
    for (let a = 0; a < 12; a++) for (let b = 0; b < 12; b++) add(6 * Math.floor(a / 3) + 2 + a % 3, 6 * Math.floor(b / 3) + 2 + b % 3, kb[a * 12 + b]);
    // drilling: γ (θz − ω)², ω = ½(∂v/∂x − ∂u/∂y), 2 x 2 Gauss points
    const gam = drill * E / (2 * (1 + nu)) * t;
    for (const [xi, eta] of GAUSS4) {
      const s = quadShape(xy, xi, eta), row = new Float64Array(24);
      for (let i = 0; i < 4; i++) { row[6 * i] = 0.5 * s.Ny[i]; row[6 * i + 1] = -0.5 * s.Nx[i]; row[6 * i + 5] = s.N[i]; }
      for (let a = 0; a < 24; a++) {
        if (row[a] === 0) continue;
        for (let b = 0; b < 24; b++) Kl[a * 24 + b] += gam * s.det * row[a] * row[b];
      }
    }
    // global: K = Tᵀ Kl T, T = blockdiag(R) for the 8 vector triplets
    const K = new Float64Array(576);
    const tmp = new Float64Array(576);
    for (let r = 0; r < 24; r++)
      for (let g = 0; g < 24; g++) {
        const bg = g - g % 3, kk = g % 3;
        let s = 0;
        for (let l = 0; l < 3; l++) s += Kl[r * 24 + bg + l] * R[l][kk];
        tmp[r * 24 + g] = s;
      }
    for (let g = 0; g < 24; g++) {
      const bg = g - g % 3, kk = g % 3;
      for (let h = 0; h < 24; h++) {
        let s = 0;
        for (let l = 0; l < 3; l++) s += R[l][kk] * tmp[(bg + l) * 24 + h];
        K[g * 24 + h] = s;
      }
    }
    const e = new LinearElement(K, 24);
    e.warp = warp;
    e.area = A;
    return e;
  }

  FA.elements = { beamIntegration, Transf, ElasticBeam, DispBeam, ForceBeam, ZeroLength, ZeroLengthRot,
                  LinearElement, quad2D, shell3D, quadArea, cross, unit };
})(globalThis.FA = globalThis.FA || {});
