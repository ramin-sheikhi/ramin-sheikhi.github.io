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

  /** Rotational spring on dof `dof` (0-based within a node) between nodes a and b. */
  class ZeroLengthRot {
    constructor(ndf, dof, makeMat) {
      this.ndf = ndf;
      this.n = 2 * ndf;
      this.dof = dof;
      this.mat = makeMat();
      this.noRayleigh = true;          // zeroLength elements carry no Rayleigh damping by default
      this.ug = new Float64Array(this.n);
      this.Cug = new Float64Array(this.n);
    }
    update(ug) {
      this.ug.set(ug);
      this.mat.setTrial(ug[this.ndf + this.dof] - ug[this.dof]);
      return true;
    }
    kmat(k) {
      const n = this.n, K = new Float64Array(n * n), i = this.dof, j = this.ndf + this.dof;
      K[i * n + i] = k; K[j * n + j] = k; K[i * n + j] = -k; K[j * n + i] = -k;
      return K;
    }
    K() { return this.kmat(this.mat.tan); }
    Kinit() { return this.kmat(this.mat.E0); }
    F() {
      const F = new Float64Array(this.n);
      F[this.dof] = -this.mat.sig;
      F[this.ndf + this.dof] = this.mat.sig;
      return F;
    }
    commit() { this.mat.commit(); this.Cug.set(this.ug); }
    revert() { this.mat.revert(); this.ug.set(this.Cug); }
  }

  FA.elements = { beamIntegration, Transf, ElasticBeam, DispBeam, ForceBeam, ZeroLengthRot };
})(globalThis.FA = globalThis.FA || {});
