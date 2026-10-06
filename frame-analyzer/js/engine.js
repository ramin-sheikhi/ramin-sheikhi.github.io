/* Finite element engine for 2D and 3D moment frames: model builder, gravity,
 * modal, pushover and time-history analyses (JavaScript port of engine.py).
 *
 * Units: kN, m, s, t (tonne = kN·s²/m). Section inputs in mm and MPa.
 * Axes: 2D frames lie in the X-Y plane with Y vertical. 3D frames use X and Y
 * in plan and Z vertical. Column depth h is along global X; beam depth h is
 * vertical.
 *
 * Constraints are imposed exactly by a transformation: every node dof is a
 * linear combination of independent unknowns (equations). Fixed dofs have no
 * unknowns, equalDOF / rigidDiaphragm slaves reuse the master's unknowns, and
 * the auxiliary nodes of beam-end hinges share every unknown with the column
 * node except the bending rotation.
 */
(function (FA) {
  'use strict';
  const LA = FA.linalg, MAT = FA.materials, SEC = FA.sections, EL = FA.elements;
  const G = 9.81, MPA = 1000.0;

  const tick = () => new Promise(r => setTimeout(r, 0));
  const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());

  // ========================================================================
  // Domain
  // ========================================================================
  class Domain {
    constructor(dim) {
      this.dim = dim;
      this.ndf = dim === 2 ? 3 : 6;
      this.nodes = new Map();
      this.elements = [];
      this.vars = [];                 // per unknown: {node, dof}
      this.messages = [];
      this.aM = 0; this.bK = 0;
    }

    log(msg) {
      const m = this.messages, last = m[m.length - 1];
      if (last && last.text === msg) last.count++;
      else m.push({ text: msg, count: 1 });
      if (m.length > 2000) m.splice(0, m.length - 2000);
    }

    logText(maxLines = 300) {
      const lines = this.messages.map(x => (x.count > 1 ? `${x.text}  (repeated ${x.count} times)` : x.text));
      if (lines.length > maxLines) return [`... (${lines.length - maxLines} earlier lines omitted)`].concat(lines.slice(-maxLines)).join('\n');
      return lines.join('\n');
    }

    addNode(tag, xyz) {
      const node = { tag, xyz, mass: 0, map: [] };
      this.nodes.set(tag, node);
      return node;
    }

    newVar(tag, dof) {
      this.vars.push({ node: tag, dof });
      return [[this.vars.length - 1, 1]];
    }

    /** Give every dof of node tag its own unknown unless already mapped. */
    freeNode(tag, skip = []) {
      const node = this.nodes.get(tag);
      for (let d = 0; d < this.ndf; d++) if (!skip.includes(d)) node.map[d] = this.newVar(tag, d);
      return node;
    }

    addElement(ele, a, b) {
      ele.nodes = [a, b];
      this.elements.push(ele);
      return ele;
    }

    /** Number the unknowns (RCM), build the skyline profile and work arrays. */
    finalize() {
      const nv = this.vars.length;
      const adj = Array.from({ length: nv }, () => new Set());
      const eleVars = this.elements.map(e => {
        const s = new Set();
        for (const t of e.nodes) for (const terms of this.nodes.get(t).map) for (const [v] of terms) s.add(v);
        return [...s];
      });
      for (const node of this.nodes.values()) {
        const s = new Set();
        for (const terms of node.map) for (const [v] of terms) s.add(v);
        eleVars.push([...s]);
      }
      for (const vs of eleVars) for (const a of vs) for (const b of vs) if (a !== b) adj[a].add(b);
      const perm = LA.rcm(adj.map(s => [...s]));          // perm[eq] = var
      const eqOf = new Int32Array(nv);
      perm.forEach((v, eq) => { eqOf[v] = eq; });
      for (const node of this.nodes.values())
        node.map = node.map.map(terms => terms.map(([v, c]) => [eqOf[v], c]));
      this.eqInfo = perm.map(v => this.vars[v]);
      const n = nv;
      this.neq = n;
      const first = new Int32Array(n);
      for (let j = 0; j < n; j++) first[j] = j;
      for (const vs of eleVars) {
        let mn = Infinity;
        for (const v of vs) mn = Math.min(mn, eqOf[v]);
        for (const v of vs) { const e = eqOf[v]; if (mn < first[e]) first[e] = mn; }
      }
      this.prof = new LA.Profile(first);
      // element location arrays
      for (const e of this.elements) {
        const ptr = [0], eqs = [], cs = [];
        for (const t of e.nodes) {
          const map = this.nodes.get(t).map;
          for (let d = 0; d < this.ndf; d++) {
            for (const [q, c] of map[d]) { eqs.push(q); cs.push(c); }
            ptr.push(eqs.length);
          }
        }
        e.ptr = Int32Array.from(ptr); e.eqs = Int32Array.from(eqs); e.cs = Float64Array.from(cs);
      }
      const P = this.prof;
      this.Kt = P.alloc(); this.Kf = P.alloc(); this.Kc = P.alloc(); this.Mv = P.alloc();
      this.K0 = null;
      for (const node of this.nodes.values()) {
        if (!node.mass) continue;
        for (let d = 0; d < this.dim; d++) {
          const terms = node.map[d];
          for (const [i, ci] of terms) for (const [j, cj] of terms) if (i <= j) P.add(this.Mv, i, j, node.mass * ci * cj);
        }
      }
      const z = () => new Float64Array(n);
      Object.assign(this, { U: z(), Uc: z(), V: z(), Vc: z(), A: z(), Ac: z(), Pl: z(), R: z(), tmp: z() });
      this.lamG = 0; this.ClamG = 0; this.lamL = 0; this.ClamL = 0; this.time = 0; this.Ctime = 0;
    }

    nodeDisp(tag, dof, vec) {                 // dof 1-based as OpenSees
      const terms = this.nodes.get(tag).map[dof - 1];
      vec = vec || this.U;
      let s = 0;
      for (const [q, c] of terms) s += c * vec[q];
      return s;
    }

    gather(e, vec) {
      const ug = new Float64Array(e.n), { ptr, eqs, cs } = e;
      for (let p = 0; p < e.n; p++) {
        let s = 0;
        for (let t = ptr[p]; t < ptr[p + 1]; t++) s += cs[t] * vec[eqs[t]];
        ug[p] = s;
      }
      return ug;
    }

    scatterK(arr, e, k) {
      const n = e.n, { ptr, eqs, cs } = e, diag = this.prof.diag;
      for (let p = 0; p < n; p++)
        for (let tp = ptr[p]; tp < ptr[p + 1]; tp++) {
          const I = eqs[tp], ci = cs[tp];
          for (let q = 0; q < n; q++) {
            const kv = k[p * n + q];
            if (kv === 0) continue;
            for (let tq = ptr[q]; tq < ptr[q + 1]; tq++) {
              const J = eqs[tq];
              if (I <= J) arr[diag[J] - (J - I)] += ci * cs[tq] * kv;
            }
          }
        }
    }

    assembleK(arr, initial, rayleighOnly) {
      arr.fill(0);
      for (const e of this.elements) {
        if (rayleighOnly && e.noRayleigh) continue;
        this.scatterK(arr, e, initial ? e.Kinit() : e.K());
      }
    }

    /** out = internal resisting forces. */
    formFint(out) {
      out.fill(0);
      for (const e of this.elements) {
        const F = e.F(), { ptr, eqs, cs } = e;
        for (let p = 0; p < e.n; p++) {
          const f = F[p];
          if (f === 0) continue;
          for (let t = ptr[p]; t < ptr[p + 1]; t++) out[eqs[t]] += cs[t] * f;
        }
      }
      return out;
    }

    updateElements() {
      for (const e of this.elements) {
        if (!e.update(this.gather(e, this.U), this.lamG)) {
          this.log(e.failMsg || 'WARNING element state determination failed');
          return false;
        }
      }
      return true;
    }

    commit() {
      for (const e of this.elements) e.commit();
      this.Uc.set(this.U); this.Vc.set(this.V); this.Ac.set(this.A);
      this.ClamG = this.lamG; this.ClamL = this.lamL; this.Ctime = this.time;
    }

    revert() {
      for (const e of this.elements) e.revert();
      this.U.set(this.Uc); this.V.set(this.Vc); this.A.set(this.Ac);
      this.lamG = this.ClamG; this.lamL = this.ClamL; this.time = this.Ctime;
    }

    nodalLoad(tag, dof, value, out) {
      for (const [q, c] of this.nodes.get(tag).map[dof - 1]) out[q] += c * value;
    }

    // ----------------------------------------------------------------------
    // Integrators
    // ----------------------------------------------------------------------
    setLoadControl(dlam) { this.mode = 'load'; this.dlam = dlam; }

    setDispControl(tag, dof, dU) {
      this.mode = 'disp';
      this.ctrl = this.nodes.get(tag).map[dof - 1];
      this.dU = dU;
      this.Kfac = this.Kf;
      if (!this.ctrl.length) throw new Error('DisplacementControl: control dof is fixed');
    }

    setNewmark(gamma, beta, accel) {
      this.mode = 'newmark';
      this.gamma = gamma; this.beta = beta;
      this.accel = accel;          // {values, dt, factor, Mr}
    }

    ctrlValue(vec) {
      let s = 0;
      for (const [q, c] of this.ctrl) s += c * vec[q];
      return s;
    }

    startStep(dt) {
      if (this.mode === 'load') {
        this.lamG += this.dlam;
        return this.updateElements();
      }
      if (this.mode === 'disp') {
        this.target = this.ctrlValue(this.U) + this.dU;
        return true;
      }
      // Newmark predictor (displacement form, as OpenSees)
      const { gamma, beta } = this;
      this.dt = dt;
      const a1 = 1 - gamma / beta, a2 = dt * (1 - 0.5 * gamma / beta);
      const a3 = -1 / (beta * dt), a4 = 1 - 0.5 / beta;
      const { V, A, Vc, Ac } = this;
      for (let i = 0; i < this.neq; i++) {
        V[i] = a1 * Vc[i] + a2 * Ac[i];
        A[i] = a4 * Ac[i] + a3 * Vc[i];
      }
      this.c2 = gamma / (beta * dt);
      this.c3 = 1 / (beta * dt * dt);
      this.time += dt;
      this.ag = this.groundAccel(this.time);
      // committed tangent for stiffness-proportional damping
      if (this.bK) this.assembleK(this.Kc, false, true);
      return true;
    }

    groundAccel(t) {
      const { values, dt, factor } = this.accel;
      // OpenSees Path series: linear interpolation, zero after the last interval
      const x = t / dt, i = Math.floor(x);
      if (i + 1 >= values.length || i < 0) return 0;
      return factor * (values[i] + (x - i) * (values[i + 1] - values[i]));
    }

    /** Effective tangent into this.Kf (not yet factored). */
    formTangent(initial) {
      const P = this.prof;
      let base;
      if (initial) {
        if (!this.K0) { this.K0 = P.alloc(); this.assembleK(this.K0, true); }
        base = this.K0;
      } else {
        this.assembleK(this.Kt, false);
        base = this.Kt;
      }
      const Kf = this.Kf;
      Kf.set(base);
      if (this.mode === 'newmark') {
        const cm = this.c3 + this.c2 * this.aM, ck = this.c2 * this.bK, M = this.Mv, Kc = this.Kc;
        for (let i = 0; i < Kf.length; i++) Kf[i] += cm * M[i] + (ck ? ck * Kc[i] : 0);
      }
    }

    factorTangent(initial) {
      this.formTangent(initial);
      try {
        this.prof.factor(this.Kf);
        return true;
      } catch (err) {
        if (err instanceof LA.SingularError) {
          const info = this.eqInfo[err.eq];
          this.log(`WARNING singular stiffness matrix - zero pivot at equation ${err.eq} ` +
            `(node ${info.node}, dof ${info.dof + 1}); LinSOE::solve failed`);
          return false;
        }
        throw err;
      }
    }

    residual(R) {
      this.formFint(R);
      const n = this.neq;
      if (this.mode === 'load') {
        for (let i = 0; i < n; i++) R[i] = -R[i];
      } else if (this.mode === 'disp') {
        const lam = this.lamL, Pl = this.Pl;
        for (let i = 0; i < n; i++) R[i] = lam * Pl[i] - R[i];
      } else {
        const P = this.prof, tmp = this.tmp, ag = this.ag, Mr = this.accel.Mr;
        for (let i = 0; i < n; i++) R[i] = -ag * Mr[i] - R[i];
        P.matvec(this.Mv, this.A, tmp);
        for (let i = 0; i < n; i++) R[i] -= tmp[i];
        if (this.aM) {
          P.matvec(this.Mv, this.V, tmp);
          for (let i = 0; i < n; i++) R[i] -= this.aM * tmp[i];
        }
        if (this.bK) {
          P.matvec(this.Kc, this.V, tmp);
          for (let i = 0; i < n; i++) R[i] -= this.bK * tmp[i];
        }
      }
      return R;
    }

    /** Correction from the factored tangent: dU (and dλ for DisplacementControl in x[n]). */
    correction(R) {
      const n = this.neq, x = new Float64Array(n + 1);
      x.set(R);
      this.prof.solve(this.Kf, x.subarray(0, n));
      if (this.mode === 'disp') {
        const b = Float64Array.from(this.Pl);
        this.prof.solve(this.Kf, b);
        const bc = this.ctrlValue(b);
        const dl = (this.target - this.ctrlValue(this.U) - this.ctrlValue(x)) / bc;
        for (let i = 0; i < n; i++) x[i] += dl * b[i];
        x[n] = dl;
      }
      return x;
    }

    apply(x, s) {
      const n = this.neq, U = this.U;
      for (let i = 0; i < n; i++) U[i] += s * x[i];
      if (this.mode === 'disp') this.lamL += s * x[n];
      else if (this.mode === 'newmark') {
        const a = s * this.c2, b = s * this.c3, V = this.V, A = this.A;
        for (let i = 0; i < n; i++) { V[i] += a * x[i]; A[i] += b * x[i]; }
      }
      return this.updateElements();
    }

    /**
     * One analysis step with the given solution algorithm. Returns true on
     * convergence (state committed) or false (state reverted).
     */
    step(alg, test, dt) {
      const { tol, maxIter } = test;
      const n = this.neq;
      const ok = this.solve(alg, tol, maxIter, dt);
      if (ok) { this.commit(); return true; }
      this.revert();
      return false;
    }

    solve(alg, tol, maxIter, dt) {
      if (!this.startStep(dt)) return false;
      const n = this.neq, R = this.R;
      const initial = alg === 'NewtonInitial';
      const everyIter = alg === 'Newton' || alg === 'NewtonLineSearch';
      const Vk = [], AVk = [];
      let factored = false;
      for (let it = 0; it < maxIter; it++) {
        if (everyIter || !factored) {
          if (!this.factorTangent(initial)) return false;
          factored = true;
        }
        this.residual(R);
        let x = this.correction(R);
        if (alg === 'KrylovNewton') {
          if (Vk.length) {
            const prev = Vk[Vk.length - 1];
            const av = new Float64Array(n + 1);
            for (let i = 0; i <= n; i++) av[i] = prev[i] - x[i];
            AVk.push(av);
            const c = LA.lstsq(AVk, x);
            const y = Float64Array.from(x);
            for (let k = 0; k < c.length; k++)
              for (let i = 0; i <= n; i++) y[i] += c[k] * (Vk[k][i] - AVk[k][i]);
            x = y;
            if (AVk.length >= 3) { Vk.length = 0; AVk.length = 0; }
          }
          Vk.push(x);
        }
        let s = 1;
        if (alg === 'NewtonLineSearch') {
          const r = this.lineSearch(x, R);
          if (r === null) return false;
          s = r;
        } else if (!this.apply(x, 1)) return false;
        const nrm = s * LA.norm(x, n);
        if (!isFinite(nrm)) return false;
        if (nrm <= tol) return true;
      }
      this.log(`WARNING ${alg}::solveCurrentStep() - failed to converge after ${maxIter} iterations ` +
        `(NormDispIncr tol ${tol}) at ${this.where()}`);
      return false;
    }

    where() {
      if (this.mode === 'newmark') return `time ${this.time.toFixed(4)}`;
      if (this.mode === 'disp') return `control displacement ${this.target.toExponential(4)}`;
      return `load factor ${this.lamG.toFixed(3)}`;
    }

    /** Interpolated line search (ratio 0.8). Applies the update; returns η. */
    lineSearch(x, R0) {
      const n = this.neq;
      const s0 = LA.dot(x.subarray(0, n), R0);
      if (!this.apply(x, 1)) return null;
      const R = this.residual(new Float64Array(n));
      let sJ = LA.dot(x.subarray(0, n), R);
      if (!(Math.abs(s0) > 0) || Math.abs(sJ) <= 0.8 * Math.abs(s0)) return 1;
      let etaPrev = 0, sU = s0, eta = 1;
      for (let i = 0; i < 10; i++) {
        let etaNew = eta - sJ * (eta - etaPrev) / (sJ - sU);
        if (!isFinite(etaNew)) break;
        etaNew = Math.min(10, Math.max(0.1, etaNew));
        if (!this.apply(x, etaNew - eta)) return null;
        etaPrev = eta; sU = sJ; eta = etaNew;
        this.residual(R);
        sJ = LA.dot(x.subarray(0, n), R);
        if (Math.abs(sJ) <= 0.8 * Math.abs(s0)) break;
      }
      return eta;
    }

    // ----------------------------------------------------------------------
    // Eigen analysis: K φ = λ M φ with the current tangent K
    // ----------------------------------------------------------------------
    eigen(nModes) {
      const n = this.neq, P = this.prof;
      this.assembleK(this.Kt, false);
      this.Kf.set(this.Kt);
      try { P.factor(this.Kf); } catch (err) {
        if (err instanceof LA.SingularError) {
          this.log(`WARNING eigen analysis: singular stiffness matrix at equation ${err.eq}`);
          throw new Error('eigen analysis failed: the stiffness matrix is singular');
        }
        throw err;
      }
      const massEqs = [];
      for (let j = 0; j < n; j++) if (this.Mv[P.diag[j]] > 0) massEqs.push(j);
      const nm = massEqs.length;
      if (nModes > nm) {
        this.log(`WARNING eigen analysis: NCV must be greater than the number of modes - ` +
          `${nModes} modes requested, the model has ${nm} dynamic dofs`);
        throw new Error(`eigen analysis failed: ${nModes} modes requested but the model has only ${nm} mass dofs`);
      }
      const res = (nm <= 1200 || nModes > nm / 3) ? this.eigenDense(massEqs, nModes) : this.eigenSubspace(massEqs, nModes);
      // sign convention: largest component positive
      for (const v of res.vectors) {
        let big = 0;
        for (let i = 0; i < n; i++) if (Math.abs(v[i]) > Math.abs(big)) big = v[i];
        if (big < 0) for (let i = 0; i < n; i++) v[i] = -v[i];
      }
      return res;
    }

    eigenDense(massEqs, nModes) {
      const n = this.neq, P = this.prof, nm = massEqs.length;
      // flexibility at the mass dofs: F = (K⁻¹)_mm = (condensed K)⁻¹
      const X = [];
      for (let j = 0; j < nm; j++) {
        const x = new Float64Array(n);
        x[massEqs[j]] = 1;
        P.solve(this.Kf, x);
        X.push(x);
      }
      const F = new Float64Array(nm * nm), Mm = new Float64Array(nm * nm);
      for (let i = 0; i < nm; i++)
        for (let j = 0; j < nm; j++) {
          F[i * nm + j] = 0.5 * (X[j][massEqs[i]] + X[i][massEqs[j]]);
          Mm[i * nm + j] = P.get(this.Mv, massEqs[i], massEqs[j]);
        }
      const L = LA.cholesky(Mm, nm);
      if (!L) throw new Error('eigen analysis failed: mass matrix is not positive definite');
      // A = Lᵀ F L
      const FL = new Float64Array(nm * nm);
      for (let i = 0; i < nm; i++)
        for (let j = 0; j < nm; j++) {
          let s = 0;
          for (let k = j; k < nm; k++) s += F[i * nm + k] * L[k * nm + j];
          FL[i * nm + j] = s;
        }
      const A = [];
      for (let i = 0; i < nm; i++) {
        const row = new Float64Array(nm);
        for (let j = 0; j < nm; j++) {
          let s = 0;
          for (let k = i; k < nm; k++) s += L[k * nm + i] * FL[k * nm + j];
          row[j] = s;
        }
        A.push(row);
      }
      for (let i = 0; i < nm; i++) for (let j = 0; j < i; j++) A[i][j] = A[j][i] = 0.5 * (A[i][j] + A[j][i]);
      const { values, vectors } = LA.symEig(A);
      const lams = [], vecs = [];
      for (let m = 0; m < nModes; m++) {
        const col = nm - 1 - m;
        const mu = values[col];
        // φ = L⁻ᵀ ψ
        const phi = new Float64Array(nm);
        for (let i = nm - 1; i >= 0; i--) {
          let s = vectors[i][col];
          for (let k = i + 1; k < nm; k++) s -= L[k * nm + i] * phi[k];
          phi[i] = s / L[i * nm + i];
        }
        const lam = 1 / mu;
        // full vector: x = λ K⁻¹ M φ = λ Σ_j X_j (Mm φ)_j
        const mphi = new Float64Array(nm);
        for (let i = 0; i < nm; i++) {
          let s = 0;
          for (let j = 0; j < nm; j++) s += Mm[i * nm + j] * phi[j];
          mphi[i] = s;
        }
        const x = new Float64Array(n);
        for (let j = 0; j < nm; j++) {
          const c = lam * mphi[j];
          if (c === 0) continue;
          const Xj = X[j];
          for (let i = 0; i < n; i++) x[i] += c * Xj[i];
        }
        massEqs.forEach((q, i) => { x[q] = phi[i]; });
        lams.push(lam);
        vecs.push(x);
      }
      return { values: lams, vectors: vecs };
    }

    eigenSubspace(massEqs, p) {
      const n = this.neq, P = this.prof;
      const q = Math.min(2 * p, p + 8, massEqs.length);
      const Md = massEqs.map(e => this.Mv[P.diag[e]]);
      const Kd = massEqs.map(e => this.Kt[P.diag[e]]);
      let X = [];
      const x0 = new Float64Array(n);
      massEqs.forEach((e, i) => { x0[e] = Md[i]; });
      X.push(x0);
      const order = massEqs.map((e, i) => i).sort((a, b) => Md[b] / Kd[b] - Md[a] / Kd[a]);
      let seed = 12345;
      const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
      for (let k = 0; X.length < q; k++) {
        const x = new Float64Array(n);
        if (k < q - 2) x[massEqs[order[k]]] = 1;
        else massEqs.forEach(e => { x[e] = rnd() - 0.5; });
        X.push(x);
      }
      let lamOld = null, result = null;
      for (let iter = 0; iter < 80; iter++) {
        const Y = X.map(x => P.matvec(this.Mv, x));
        const Xb = Y.map(y => P.solve(this.Kf, Float64Array.from(y)));
        const MXb = Xb.map(x => P.matvec(this.Mv, x));
        const Kr = [], Mr = [];
        for (let i = 0; i < q; i++) {
          Kr.push(new Float64Array(q)); Mr.push(new Float64Array(q));
        }
        for (let i = 0; i < q; i++)
          for (let j = 0; j <= i; j++) {
            Kr[i][j] = Kr[j][i] = 0.5 * (LA.dot(Xb[i], Y[j]) + LA.dot(Xb[j], Y[i]));
            Mr[i][j] = Mr[j][i] = 0.5 * (LA.dot(Xb[i], MXb[j]) + LA.dot(Xb[j], MXb[i]));
          }
        const { values, vectors } = LA.genEig(Kr, Mr);
        X = vectors.map(Q => {
          const x = new Float64Array(n);
          for (let k = 0; k < q; k++) {
            const c = Q[k];
            if (c === 0) continue;
            for (let i = 0; i < n; i++) x[i] += c * Xb[k][i];
          }
          return x;
        });
        result = { values: Array.from(values.slice(0, p)), vectors: X.slice(0, p) };
        if (lamOld) {
          let conv = true;
          for (let i = 0; i < p; i++) if (Math.abs(values[i] - lamOld[i]) > 1e-10 * Math.abs(values[i])) conv = false;
          if (conv) break;
        }
        lamOld = Array.from(values);
      }
      return result;
    }
  }

  // ========================================================================
  // Geometry helpers
  // ========================================================================
  const nodeTag = (floor, ix, iy = 0) => floor * 10000 + iy * 100 + ix + 1;

  /** [X grid, Y grid, level elevations]. The Y grid is [0] in 2D. */
  function frameGeometry(p) {
    const levels = [0.0];
    for (let i = 0; i < p.n_stories; i++) levels.push(levels[levels.length - 1] + (i === 0 ? p.h1 : p.h));
    const xs = [];
    for (let j = 0; j <= p.n_bays; j++) xs.push(j * p.bay);
    const yg = [];
    if ((p.dim || 2) === 3) for (let k = 0; k <= p.n_bays_y; k++) yg.push(k * p.bay_y);
    else yg.push(0.0);
    return [xs, yg, levels];
  }

  const elasticProps = s => {
    const sp = SEC.sectionProps(s);
    return [sp.E, sp.A, sp.Iz];
  };

  // ========================================================================
  // Model
  // ========================================================================
  /**
   * Create the frame (no loads analysed yet) and return a model-info object.
   * p: dim (2|3), n_stories, n_bays, bay (X), n_bays_y, bay_y (3D only),
   *    h1, h, rigid_diaphragm, pdelta, ele_type, n_ip, integ,
   *    col_sec, beam_sec, w_floor, w_roof (kN/m on every beam),
   *    m_floor, m_roof (t per floor),
   *    conn: {type: 'Rigid'|'Pinned'|'Semi-rigid', mat: spring material spec (kN·m, rad)}
   */
  function createModel(p) {
    const dim = p.dim || 2, ndf = dim === 2 ? 3 : 6;
    const D = new Domain(dim);
    const [xs, yg, levels] = frameGeometry(p);
    const ns = p.n_stories;
    const plan = [];
    for (let k = 0; k < yg.length; k++) for (let j = 0; j < xs.length; j++) plan.push([j, k]);
    const nn = plan.length;
    const point = (i, j, k) => (dim === 2 ? [xs[j], levels[i]] : [xs[j], yg[k], levels[i]]);

    const coords = {}, allNodes = [], planNodes = [];
    for (let i = 0; i <= ns; i++) {
      planNodes.push([]);
      for (const [j, k] of plan) {
        const n = nodeTag(i, j, k);
        D.addNode(n, point(i, j, k));
        coords[n] = point(i, j, k);
        allNodes.push(n);
        planNodes[i].push(n);
      }
    }
    for (const n of planNodes[0]) D.nodes.get(n).map = Array.from({ length: ndf }, () => []);

    // Masses (lumped, split equally between the nodes of each floor)
    const floorMass = [], nodeMass = {};
    for (let i = 1; i <= ns; i++) {
      const m = i === ns ? p.m_roof : p.m_floor;
      floorMass.push(m);
      for (const n of planNodes[i]) {
        nodeMass[n] = m / nn;
        D.nodes.get(n).mass = m / nn;
      }
    }

    // Free dofs and rigid floors; floorNodes[i] represents level i
    const floorNodes = [planNodes[0][0]];
    for (let i = 1; i <= ns; i++) {
      const nodes = planNodes[i];
      if (!p.rigid_diaphragm) {
        for (const n of nodes) D.freeNode(n);
        floorNodes.push(nodes[0]);
      } else if (dim === 2) {
        const m0 = D.freeNode(nodes[0]);
        for (const n of nodes.slice(1)) {
          const nd = D.freeNode(n, [0]);
          nd.map[0] = m0.map[0];                    // equalDOF(master, n, 1)
        }
        floorNodes.push(nodes[0]);
      } else {
        const master = 2000000 + i;
        const xm = (xs[0] + xs[xs.length - 1]) / 2, ym = (yg[0] + yg[yg.length - 1]) / 2;
        const mn = D.addNode(master, [xm, ym, levels[i]]);
        mn.map = [D.newVar(master, 0), D.newVar(master, 1), [], [], [], D.newVar(master, 5)];
        const [ux] = mn.map[0], [uy] = mn.map[1], [rz] = mn.map[5];
        for (const n of nodes) {
          const nd = D.freeNode(n, [0, 1, 5]);    // rigidDiaphragm(3, master, ...)
          const [x, y] = coords[n];
          nd.map[0] = [[ux[0], 1], [rz[0], -(y - ym)]];
          nd.map[1] = [[uy[0], 1], [rz[0], x - xm]];
          nd.map[5] = [[rz[0], 1]];
        }
        coords[master] = [xm, ym, levels[i]];
        floorNodes.push(master);
      }
    }
    let nDyn;
    if (dim === 2) nDyn = p.rigid_diaphragm ? ns * (1 + nn) : 2 * ns * nn;
    else nDyn = p.rigid_diaphragm ? ns * (3 + nn) : 3 * ns * nn;

    // Members
    const colPD = !!p.pdelta;
    const vec = { col: [0, 1, 0], bx: [0, -1, 0], by: [1, 0, 0] };
    const transf = (a, b, kind) => new EL.Transf(dim, coords[a], coords[b], dim === 3 ? vec[kind] : null,
      kind === 'col' && colPD);

    function member(spec) {
      if (p.ele_type === 'elasticBeamColumn') {
        const sp = SEC.sectionProps(spec);
        return (tr, wy) => new EL.ElasticBeam(dim, tr, sp, wy);
      }
      const mk = SEC.sectionMaker(spec, dim);
      const rule = EL.beamIntegration(p.integ, p.n_ip);
      const Cls = p.ele_type === 'dispBeamColumn' ? EL.DispBeam : EL.ForceBeam;
      return (tr, wy) => new Cls(dim, tr, mk, rule, wy);
    }

    const makeCol = member(p.col_sec);
    const makeBeam = member(p.beam_sec);
    let eleTag = 0;
    const columns = [];
    for (let i = 0; i < ns; i++) {
      planNodes[i].forEach((a, idx) => {
        const b = planNodes[i + 1][idx];
        const t = ++eleTag;
        const e = D.addElement(makeCol(transf(a, b, 'col'), 0), a, b);
        e.tag = t;
        columns.push([t, a, b]);
      });
    }

    const conn = p.conn || { type: 'Rigid' };
    const springs = [];
    const springMat = conn.type === 'Semi-rigid' ? MAT.maker(conn.mat, 1.0) : null;
    let aux = 1000000;

    function beamEnd(n, bdir) {
      if (conn.type === 'Rigid') return n;
      const bn = ++aux;
      D.addNode(bn, coords[n]);
      coords[bn] = coords[n];
      allNodes.push(bn);
      const bend = dim === 2 ? 2 : (bdir === 1 ? 4 : 3);
      const src = D.nodes.get(n).map;
      D.nodes.get(bn).map = src.map((terms, d) => (d === bend ? D.newVar(bn, d) : terms));
      const t = ++eleTag;
      if (springMat) {
        const e = D.addElement(new EL.ZeroLengthRot(ndf, bend, springMat), n, bn);
        e.tag = t;
      }
      springs.push([t, n, bn]);
      return bn;
    }

    const beams = [];
    for (let i = 1; i <= ns; i++) {
      const lines = [];
      for (let k = 0; k < yg.length; k++)
        for (let j = 0; j < xs.length - 1; j++) lines.push([1, nodeTag(i, j, k), nodeTag(i, j + 1, k)]);
      if (dim === 3)
        for (let j = 0; j < xs.length; j++)
          for (let k = 0; k < yg.length - 1; k++) lines.push([2, nodeTag(i, j, k), nodeTag(i, j, k + 1)]);
      const w = i === ns ? p.w_roof : p.w_floor;
      for (const [bdir, na, nb] of lines) {
        const t = ++eleTag;
        const a = beamEnd(na, bdir), b = beamEnd(nb, bdir);
        const e = D.addElement(makeBeam(transf(a, b, bdir === 1 ? 'bx' : 'by'), -w), a, b);
        e.tag = t;
        beams.push([t, a, b, i, bdir]);
      }
    }

    D.finalize();
    const baseSet = new Set(planNodes[0]);
    D.baseEls = D.elements.filter(e => baseSet.has(e.nodes[0]));

    const info = {
      dim, xs, yg, levels, ns, columns, beams, springs, conn: conn.type, floor_mass: floorMass,
      node_mass: nodeMass, base_nodes: planNodes[0], plan_nodes: planNodes, floor_nodes: floorNodes,
      ctrl_node: floorNodes[floorNodes.length - 1], coords, all_nodes: allNodes,
      max_modes: Math.max(1, nDyn), n_eq: D.neq, domain: D,
    };
    return info;
  }

  /** Build the frame and apply gravity (synchronous). */
  function buildModel(p) {
    const info = createModel(p);
    for (let k = 0; k < 10; k++) gravityStep(info, k);
    endGravity(info);
    return info;
  }

  /**
   * Same as buildModel, yielding to the browser between gravity steps.
   * hooks.cancelled() -> true aborts and resolves to null.
   */
  async function buildModelAsync(p, hooks = {}) {
    const info = createModel(p);
    for (let k = 0; k < 10; k++) {
      await tick();
      if (hooks.cancelled && hooks.cancelled()) return null;
      gravityStep(info, k);
    }
    endGravity(info);
    return info;
  }

  const ALGORITHMS = ['Newton', 'ModifiedNewton', 'NewtonLineSearch', 'KrylovNewton', 'NewtonInitial'];

  /** Gravity: LoadControl 0.1 x 10 steps, Newton, NormDispIncr 1e-8. */
  function gravityStep(info, k) {
    const D = info.domain;
    if (k === 0) D.setLoadControl(0.1);
    if (!D.step('Newton', { tol: 1e-8, maxIter: 50 })) {
      D.log(`WARNING Analysis::analyze() - gravity step ${k + 1} of 10 failed to converge`);
      const err = new Error('Gravity analysis did not converge - check sections/loads.');
      err.log = D.logText();
      throw err;
    }
  }

  function endGravity(info) {               // loadConst -time 0
    info.domain.time = 0;
    info.domain.Ctime = 0;
  }

  /** Translational displacements of every drawn node. */
  function currentDisps(info) {
    const D = info.domain, out = {};
    for (const n of info.all_nodes) {
      const u = [];
      for (let d = 1; d <= info.dim; d++) u.push(D.nodeDisp(n, d));
      out[n] = u;
    }
    return out;
  }

  const lateralDirs = info => (info.dim === 2 ? [1] : [1, 2]);

  function baseShear(info, d = 1) {
    let s = 0;
    for (const e of info.domain.baseEls) s += e.F()[d - 1];
    return -s;
  }

  const floorDisp = (info, d = 1) => info.floor_nodes.map(n => info.domain.nodeDisp(n, d));

  function tryStep(D, test, dt) {
    for (const alg of ALGORITHMS) if (D.step(alg, test, dt)) return true;
    return false;
  }

  // ========================================================================
  // Modal analysis
  // ========================================================================
  function eigenvector(info, vec, n, d) {
    return info.domain.nodeDisp(n, d, vec);
  }

  function runModal(info, nModes) {
    nModes = Math.max(1, Math.min(nModes, info.max_modes));
    const D = info.domain;
    const { values, vectors } = D.eigen(nModes);
    const dim = info.dim;
    const nodes = Object.keys(info.node_mass).map(Number);
    const m = nodes.map(n => info.node_mass[n]);
    const Mtot = m.reduce((a, b) => a + b, 0);
    const modes = [];
    for (let k = 0; k < nModes; k++) {
      const omega = Math.sqrt(Math.abs(values[k]));
      const T = 2 * Math.PI / omega;
      const v = vectors[k];
      const ph = nodes.map(n => { const r = []; for (let d = 1; d <= dim; d++) r.push(eigenvector(info, v, n, d)); return r; });
      let mn = 0;
      ph.forEach((r, i) => r.forEach(x => { mn += m[i] * x * x; }));
      const shape = {};
      for (const n of info.all_nodes) {
        const r = [];
        for (let d = 1; d <= dim; d++) r.push(eigenvector(info, v, n, d));
        shape[n] = r;
      }
      const mode = { mode: k + 1, T, f: 1 / T, omega, lambda: values[k], shape };
      for (const d of lateralDirs(info)) {
        let L = 0;
        ph.forEach((r, i) => { L += m[i] * r[d - 1]; });
        const ax = 'xy'[d - 1];
        mode['gamma_' + ax] = L / mn;
        mode['ratio_' + ax] = L * L / mn / Mtot;
      }
      modes.push(mode);
    }
    for (const d of lateralDirs(info)) {
      const ax = 'xy'[d - 1];
      let c = 0;
      for (const md of modes) { c += md['ratio_' + ax]; md['cum_' + ax] = c; }
    }
    return modes;
  }

  // ========================================================================
  // Pushover
  // ========================================================================
  function firstModeShape(info, d) {
    const D = info.domain;
    const n = Math.min(6, info.max_modes);
    const { vectors } = D.eigen(n);
    const nodes = Object.keys(info.node_mass).map(Number);
    const m = nodes.map(x => info.node_mass[x]);
    let best = 0, bestR = -1;
    vectors.forEach((v, k) => {
      let L = 0, mn = 0;
      nodes.forEach((x, i) => {
        for (let j = 1; j <= info.dim; j++) { const e = D.nodeDisp(x, j, v); mn += m[i] * e * e; }
        L += m[i] * D.nodeDisp(x, d, v);
      });
      const r = L * L / mn;
      if (r > bestR) { best = k; bestR = r; }
    });
    const phi = info.floor_nodes.slice(1).map(x => D.nodeDisp(x, d, vectors[best]));
    const sg = Math.sign(phi[phi.length - 1]) || 1;
    return phi.map(v => v * sg);
  }

  async function runPushover(info, pattern, targetDrift, nSteps, direction = 1, hooks = {}) {
    const d = direction, D = info.domain, levels = info.levels;
    const m = info.floor_mass, h = levels.slice(1);
    let w;
    if (pattern === 'Triangular (inverted)') w = m.map((mi, i) => mi * h[i]);
    else if (pattern === 'Uniform') w = m.slice();
    else { const phi = firstModeShape(info, d); w = m.map((mi, i) => mi * phi[i]); }
    const ws = w.reduce((a, b) => a + b, 0);
    w = w.map(x => x / ws);

    D.Pl.fill(0);
    for (let i = 1; i <= info.ns; i++) {
      const nodes = info.plan_nodes[i];
      for (const n of nodes) D.nodalLoad(n, d, w[i - 1] / nodes.length, D.Pl);
    }
    const H = levels[levels.length - 1];
    const target = targetDrift / 100 * H;
    const dU = target / nSteps;
    const ctrl = info.ctrl_node;
    const test = { tol: 1e-6, maxIter: 100 };
    D.setDispControl(ctrl, d, dU);

    const roof = [0], V = [0], floors = [floorDisp(info, d)];
    let status = 'Reached target displacement';
    let last = now();
    while (D.nodeDisp(ctrl, d) < target - 1e-12) {
      let ok = tryStep(D, test);
      if (!ok) {
        D.dU = dU / 10;
        ok = true;
        for (let s = 0; s < 10 && ok; s++) ok = tryStep(D, test);
        D.dU = dU;
      }
      if (!ok) {
        status = `Stopped at roof drift ${(100 * D.nodeDisp(ctrl, d) / H).toFixed(2)}% (no convergence)`;
        break;
      }
      roof.push(D.nodeDisp(ctrl, d));
      V.push(baseShear(info, d));
      floors.push(floorDisp(info, d));
      if (now() - last > 40) {
        if (hooks.progress) hooks.progress(Math.min(1, roof[roof.length - 1] / target));
        await tick();
        last = now();
        if (hooks.cancelled && hooks.cancelled()) {
          status = `Stopped by user at roof drift ${(100 * D.nodeDisp(ctrl, d) / H).toFixed(2)}%`;
          break;
        }
      }
    }
    const drifts = floors.map(f => f.slice(1).map((x, i) => (x - f[i]) / (levels[i + 1] - levels[i])));
    return { roof, V, floors, drifts, H, W: m.reduce((a, b) => a + b, 0) * G, status, final: currentDisps(info) };
  }

  // ========================================================================
  // Time history
  // ========================================================================
  async function runTimeHistory(info, accG, dtRec, scale, zeta, dtAn, extraTime = 0, dampModes = [1, 3],
    direction = 1, hooks = {}) {
    const d = direction, D = info.domain;
    const nmax = info.max_modes;
    const [mi, mj] = dampModes.slice().sort((a, b) => a - b).map(x => Math.max(1, Math.min(x, nmax)));
    const lam = D.eigen(Math.max(mi, mj)).values;
    const wi = Math.sqrt(lam[mi - 1]), wj = Math.sqrt(lam[mj - 1]);
    let aM, bK;
    if (Math.abs(wj - wi) < 1e-9) { aM = zeta * wi; bK = zeta / wi; }
    else { aM = 2 * zeta * wi * wj / (wi + wj); bK = 2 * zeta / (wi + wj); }
    D.aM = aM; D.bK = bK;          // rayleigh(aM, 0, 0, bK): committed stiffness

    // influence vector of the excitation direction and M·r
    const r = new Float64Array(D.neq);
    D.eqInfo.forEach((v, q) => { if (v.dof === d - 1) r[q] = 1; });
    const Mr = D.prof.matvec(D.Mv, r);
    D.setNewmark(0.5, 0.25, { values: accG, dt: dtRec, factor: scale * G, Mr });
    const test = { tol: 1e-8, maxIter: 50 };

    const duration = accG.length * dtRec + extraTime;
    const n = Math.round(duration / dtAn);
    const t = [0], roof = [0], V = [0], floors = [floorDisp(info, d)];
    let status = 'Completed';
    let last = now();
    for (let k = 0; k < n; k++) {
      let ok = tryStep(D, test, dtAn);
      if (!ok) {
        ok = true;
        for (let s = 0; s < 10 && ok; s++) ok = tryStep(D, test, dtAn / 10);
      }
      if (!ok) {
        status = `Stopped at t = ${D.time.toFixed(3)} s (no convergence)`;
        break;
      }
      t.push(D.time);
      roof.push(D.nodeDisp(info.ctrl_node, d));
      V.push(baseShear(info, d));
      floors.push(floorDisp(info, d));
      if (now() - last > 40) {
        if (hooks.progress) hooks.progress((k + 1) / n);
        await tick();
        last = now();
        if (hooks.cancelled && hooks.cancelled()) {
          status = `Stopped by user at t = ${D.time.toFixed(3)} s`;
          break;
        }
      }
    }
    const levels = info.levels;
    const drifts = floors.map(f => f.slice(1).map((x, i) => (x - f[i]) / (levels[i + 1] - levels[i])));
    return { t, roof, V, floors, drifts, status, rayleigh: [aM, bK],
             periods: [2 * Math.PI / wi, 2 * Math.PI / wj], W: info.floor_mass.reduce((a, b) => a + b, 0) * G };
  }

  FA.engine = {
    G, MPA, Domain, nodeTag, frameGeometry, elasticProps, createModel, buildModel, buildModelAsync, currentDisps, lateralDirs,
    baseShear, floorDisp, runModal, runPushover, runTimeHistory, ALGORITHMS,
    sectionProps: SEC.sectionProps, rcBarLayout: SEC.rcBarLayout, rcSteelSummary: SEC.rcSteelSummary,
  };
})(globalThis.FA = globalThis.FA || {});
