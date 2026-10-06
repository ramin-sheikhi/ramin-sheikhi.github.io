/* Uniaxial materials and the catalogue offered in the interface.
 *
 * Part 1 - material models. Each object has
 *     setTrial(strain)   computes the trial stress .sig and tangent .tan
 *                        from the last committed state
 *     commit()           accepts the trial state
 *     revert()           returns to the last committed state
 * Signs follow OpenSees: tension positive, concrete parameters negative.
 *
 * Part 2 - catalogue (port of materials.py). Values are entered in "input
 * units" and the force-like ones are multiplied by a scale factor:
 *     kind "f"    force / stress       MPa in sections, kN·m in joint springs
 *     kind "k"    stiffness / modulus  MPa in sections, kN·m/rad in joint springs
 *     kind "e"    strain / rotation    -, rad
 *     kind "r"    dimensionless
 *     kind "env"  HystereticSM envelope {pos: [[F, e], ...], neg: [...] | null}
 *     kind "text" free text (Custom material)
 */
(function (FA) {
  'use strict';

  const EPS = 2.220446049250313e-16;

  // ========================================================================
  // Part 1 - material models
  // ========================================================================
  class Elastic {
    constructor(E) {
      this.E = E;
      this.E0 = E;
      this.eps = 0; this.sig = 0; this.tan = E;
      this.Ceps = 0;
    }
    setTrial(e) { this.eps = e; this.sig = this.E * e; this.tan = this.E; }
    commit() { this.Ceps = this.eps; }
    revert() { this.setTrial(this.Ceps); }
  }

  /** Elastic, no tension (ENT). */
  class ENT {
    constructor(E) {
      this.E = E;
      this.E0 = E;
      this.eps = 0; this.sig = 0; this.tan = E;
      this.Ceps = 0;
    }
    setTrial(e) {
      this.eps = e;
      if (e < 0) { this.sig = this.E * e; this.tan = this.E; }
      else { this.sig = 0; this.tan = 0; }
    }
    commit() { this.Ceps = this.eps; }
    revert() { this.setTrial(this.Ceps); }
  }

  /** Elastic-perfectly plastic (OpenSees ElasticPP). */
  class ElasticPP {
    constructor(E, epsyP, epsyN, eps0) {
      this.E = E;
      this.E0 = E;
      this.fyp = E * epsyP;
      this.fyn = E * (epsyN === undefined ? -epsyP : epsyN);
      this.ez = eps0 || 0;
      this.Cep = 0; this.Ceps = 0;
      this.eps = 0; this.sig = 0; this.tan = E;
    }
    setTrial(e) {
      this.eps = e;
      const st = this.E * (e - this.ez - this.Cep);
      const f = st >= 0 ? st - this.fyp : -st + this.fyn;
      if (f <= -this.E * EPS) { this.sig = st; this.tan = this.E; }
      else { this.sig = st > 0 ? this.fyp : this.fyn; this.tan = 0; }
    }
    commit() {
      if (this.sig >= this.fyp) this.Cep = this.eps - this.ez - this.fyp / this.E;
      else if (this.sig <= this.fyn) this.Cep = this.eps - this.ez - this.fyn / this.E;
      this.Ceps = this.eps;
    }
    revert() { this.setTrial(this.Ceps); }
  }

  /** Bilinear steel with kinematic (and optional isotropic) hardening (Steel01). */
  function cpSteel01(d, s) { d.min = s.min; d.max = s.max; d.sP = s.sP; d.sN = s.sN; d.load = s.load; d.eps = s.eps; d.sig = s.sig; d.tan = s.tan; return d; }

  class Steel01 {
    constructor(fy, E0, b, a1 = 0, a2 = 55, a3 = 0, a4 = 55) {
      Object.assign(this, { fy, E0, b, a1, a2, a3, a4 });
      this.C = { min: 0, max: 0, sP: 1, sN: 1, load: 0, eps: 0, sig: 0, tan: E0 };
      this.T = cpSteel01({}, this.C);
      this.sig = 0; this.tan = E0; this.eps = 0;
    }
    setTrial(strain) {
      const C = this.C, T = this.T;
      cpSteel01(T, C);
      const dS = strain - C.eps;
      if (Math.abs(dS) > EPS) {
        T.eps = strain;
        const { fy, E0, b } = this;
        const fyOneMinusB = fy * (1 - b), Esh = b * E0, epsy = fy / E0;
        if (T.load === 0 && dS !== 0) T.load = dS > 0 ? 1 : -1;
        if (T.load === 1 && dS < 0) {
          T.load = -1;
          if (C.eps > T.max) T.max = C.eps;
          T.sN = 1 + this.a1 * Math.pow((T.max - T.min) / (2 * this.a2 * epsy), 0.8);
        }
        if (T.load === -1 && dS > 0) {
          T.load = 1;
          if (C.eps < T.min) T.min = C.eps;
          T.sP = 1 + this.a3 * Math.pow((T.max - T.min) / (2 * this.a4 * epsy), 0.8);
        }
        const c1 = Esh * strain, c2 = T.sN * fyOneMinusB, c3 = T.sP * fyOneMinusB;
        const c = C.sig + E0 * dS;
        const c1c3 = c1 + c3;
        T.sig = c1c3 < c ? c1c3 : c;
        const c1c2 = c1 - c2;
        if (c1c2 > T.sig) T.sig = c1c2;
        T.tan = Math.abs(T.sig - c) < EPS ? E0 : Esh;
      }
      this.eps = T.eps; this.sig = T.sig; this.tan = T.tan;
    }
    commit() { cpSteel01(this.C, this.T); }
    revert() { cpSteel01(this.T, this.C); this.eps = this.C.eps; this.sig = this.C.sig; this.tan = this.C.tan; }
  }

  /**
   * Giuffré-Menegotto-Pinto steel (Steel02). With different yield strengths
   * and hardening ratios in tension and compression it is also used for
   * SteelMPF.
   */
  function cpSteel02(d, s) { d.epsmax = s.epsmax; d.epsmin = s.epsmin; d.epspl = s.epspl; d.epss0 = s.epss0; d.sigs0 = s.sigs0; d.epsr = s.epsr; d.sigr = s.sigr; d.kon = s.kon; d.eps = s.eps; d.sig = s.sig; d.e = s.e; return d; }

  class Steel02 {
    constructor(fyp, fyn, E0, bp, bn, R0, cR1, cR2, a1 = 0, a2 = 1, a3 = 0, a4 = 1, sigini = 0) {
      Object.assign(this, { fyp, fyn, E0, bp, bn, R0, cR1, cR2, a1, a2, a3, a4, sigini });
      this.epsyp = fyp / E0;
      this.epsyn = fyn / E0;
      this.C = { epsmax: this.epsyp, epsmin: -this.epsyn, epspl: 0, epss0: 0, sigs0: 0,
                 epsr: 0, sigr: 0, kon: 0, eps: 0, sig: 0, e: E0 };
      this.T = cpSteel02({}, this.C);
      this.sig = 0; this.tan = E0; this.eps = 0;
    }
    setTrial(trialStrain) {
      const C = this.C, T = this.T;
      cpSteel02(T, C);
      const E0 = this.E0;
      const eps = trialStrain + this.sigini / E0;
      const deps = eps - C.eps;
      T.eps = eps;
      if (T.kon === 0 || T.kon === 3) {
        if (Math.abs(deps) < 10 * EPS) {
          T.e = E0; T.sig = this.sigini; T.kon = 3;
          this.eps = trialStrain; this.sig = T.sig; this.tan = T.e;
          return;
        }
        T.epsmax = this.epsyp;
        T.epsmin = -this.epsyn;
        if (deps < 0) { T.kon = 2; T.epss0 = T.epsmin; T.sigs0 = -this.fyn; T.epspl = T.epsmin; }
        else { T.kon = 1; T.epss0 = T.epsmax; T.sigs0 = this.fyp; T.epspl = T.epsmax; }
      }
      if (T.kon === 2 && deps > 0) {
        T.kon = 1;
        T.epsr = C.eps; T.sigr = C.sig;
        if (C.eps < T.epsmin) T.epsmin = C.eps;
        const d1 = (T.epsmax - T.epsmin) / (2 * this.a4 * this.epsyp);
        const shft = 1 + this.a3 * Math.pow(d1, 0.8);
        const Esh = this.bp * E0, fy = this.fyp, ey = this.epsyp;
        T.epss0 = (fy * shft - Esh * ey * shft - T.sigr + E0 * T.epsr) / (E0 - Esh);
        T.sigs0 = fy * shft + Esh * (T.epss0 - ey * shft);
        T.epspl = T.epsmax;
      } else if (T.kon === 1 && deps < 0) {
        T.kon = 2;
        T.epsr = C.eps; T.sigr = C.sig;
        if (C.eps > T.epsmax) T.epsmax = C.eps;
        const d1 = (T.epsmax - T.epsmin) / (2 * this.a2 * this.epsyn);
        const shft = 1 + this.a1 * Math.pow(d1, 0.8);
        const Esh = this.bn * E0, fy = this.fyn, ey = this.epsyn;
        T.epss0 = (-fy * shft + Esh * ey * shft - T.sigr + E0 * T.epsr) / (E0 - Esh);
        T.sigs0 = -fy * shft + Esh * (T.epss0 + ey * shft);
        T.epspl = T.epsmin;
      }
      const b = T.kon === 1 ? this.bp : this.bn;
      const epsy = T.kon === 1 ? this.epsyp : this.epsyn;
      const xi = Math.abs((T.epspl - T.epss0) / epsy);
      const R = this.R0 * (1 - (this.cR1 * xi) / (this.cR2 + xi));
      const epsrat = (eps - T.epsr) / (T.epss0 - T.epsr);
      const dum1 = 1 + Math.pow(Math.abs(epsrat), R);
      const dum2 = Math.pow(dum1, 1 / R);
      T.sig = (b * epsrat + (1 - b) * epsrat / dum2) * (T.sigs0 - T.sigr) + T.sigr;
      T.e = (b + (1 - b) / (dum1 * dum2)) * (T.sigs0 - T.sigr) / (T.epss0 - T.epsr);
      this.eps = trialStrain; this.sig = T.sig; this.tan = T.e;
    }
    commit() { cpSteel02(this.C, this.T); }
    revert() {
      cpSteel02(this.T, this.C);
      this.eps = this.C.eps - this.sigini / this.E0; this.sig = this.C.sig; this.tan = this.C.e;
    }
  }

  /**
   * Hysteretic material with a multi-point envelope (Hysteretic: 2-3 points,
   * HystereticSM: 2-7 points). Pinching, ductility/energy damage and degraded
   * unloading stiffness follow OpenSees HystereticMaterial.
   * pos, neg: arrays of [stress, strain] with positive magnitudes.
   */
  function cpHysteretic(d, s) { d.rotMax = s.rotMax; d.rotMin = s.rotMin; d.rotPu = s.rotPu; d.rotNu = s.rotNu; d.energyD = s.energyD; d.load = s.load; d.stress = s.stress; d.strain = s.strain; d.tangent = s.tangent; return d; }

  class Hysteretic {
    constructor(pos, neg, pinchX, pinchY, damfc1, damfc2, beta) {
      Object.assign(this, { pinchX, pinchY, damfc1, damfc2, beta: beta || 0 });
      this.P = envelope(pos);
      this.N = envelope(neg);
      this.rot1p = this.P.e[0]; this.rot1n = -this.N.e[0];
      this.E1p = this.P.E[0]; this.E1n = this.N.E[0];
      this.E0 = this.E1p;
      this.energyA = this.P.area + this.N.area;
      this.C = { rotMax: 0, rotMin: 0, rotPu: 0, rotNu: 0, energyD: 0, load: 0,
                 stress: 0, strain: 0, tangent: this.E1p };
      this.T = cpHysteretic({}, this.C);
      this.sig = 0; this.tan = this.E1p; this.eps = 0;
    }
    posEnvStress(x) { return envStress(this.P, x); }
    negEnvStress(x) { return -envStress(this.N, -x); }
    posEnvTangent(x) { return envTangent(this.P, x); }
    negEnvTangent(x) { return envTangent(this.N, -x); }
    /** Deformation at which a softening positive envelope has lost all strength (+∞ if none). */
    posEnvRotlim(x) {
      const lim = envRotlim(this.P, x);
      return lim === null ? Infinity : lim;
    }
    negEnvRotlim(x) {
      const lim = envRotlim(this.N, -x);
      return lim === null ? -Infinity : -lim;
    }
    setTrial(strain) {
      const C = this.C, T = this.T;
      cpHysteretic(T, C);
      T.strain = strain;
      const dStrain = strain - C.strain;
      if (T.load === 0) T.load = dStrain < 0 ? 2 : 1;
      if (strain >= C.rotMax) {
        T.rotMax = strain;
        T.tangent = this.posEnvTangent(strain);
        T.stress = this.posEnvStress(strain);
      } else if (strain <= C.rotMin) {
        T.rotMin = strain;
        T.tangent = this.negEnvTangent(strain);
        T.stress = this.negEnvStress(strain);
      } else if (dStrain < 0) this.negativeIncrement(dStrain);
      else if (dStrain > 0) this.positiveIncrement(dStrain);
      T.energyD = C.energyD + 0.5 * (C.stress + T.stress) * dStrain;
      this.eps = strain; this.sig = T.stress; this.tan = T.tangent;
    }
    positiveIncrement(dStrain) {
      const C = this.C, T = this.T;
      let kn = Math.pow(C.rotMin / this.rot1n, this.beta);
      kn = kn < 1 ? 1 : 1 / kn;
      let kp = Math.pow(C.rotMax / this.rot1p, this.beta);
      kp = kp < 1 ? 1 : 1 / kp;
      const Eup = this.E1p, Eun = this.E1n;
      if (T.load === 2) {
        T.load = 1;
        if (C.stress <= 0) {
          T.rotNu = C.strain - C.stress / (Eun * kn);
          const energy = C.energyD - 0.5 * C.stress / (Eun * kn) * C.stress;
          let damfac = 0;
          if (C.rotMin < this.rot1n) {
            damfac = this.damfc2 * energy / this.energyA;
            damfac += this.damfc1 * (C.rotMin / this.rot1n - 1);
          }
          T.rotMax = C.rotMax * (1 + damfac);
        }
      }
      T.load = 1;
      T.rotMax = T.rotMax > this.rot1p ? T.rotMax : this.rot1p;
      const maxmom = this.posEnvStress(T.rotMax);
      const rotlim = this.negEnvRotlim(C.rotMin);
      const rotrel = rotlim > T.rotNu ? rotlim : T.rotNu;
      const rotmp2 = T.rotMax - (1 - this.pinchY) * maxmom / (Eup * kp);
      const rotch = rotrel + (rotmp2 - rotrel) * this.pinchX;
      const s = T.strain;
      if (s < T.rotNu) {
        T.tangent = Eun * kn;
        T.stress = C.stress + T.tangent * dStrain;
        if (T.stress >= 0) { T.stress = 0; T.tangent = Eun * 1e-9; }
      } else if (s >= T.rotNu && s < rotch) {
        if (s <= rotrel) { T.stress = 0; T.tangent = Eup * 1e-9; }
        else {
          T.tangent = maxmom * this.pinchY / (rotch - rotrel);
          const t1 = C.stress + Eup * kp * dStrain;
          const t2 = (s - rotrel) * T.tangent;
          if (t1 < t2) { T.stress = t1; T.tangent = Eup * kp; } else T.stress = t2;
        }
      } else {
        T.tangent = (1 - this.pinchY) * maxmom / (T.rotMax - rotch);
        const t1 = C.stress + Eup * kp * dStrain;
        const t2 = this.pinchY * maxmom + (s - rotch) * T.tangent;
        if (t1 < t2) { T.stress = t1; T.tangent = Eup * kp; } else T.stress = t2;
      }
    }
    negativeIncrement(dStrain) {
      const C = this.C, T = this.T;
      let kn = Math.pow(C.rotMin / this.rot1n, this.beta);
      kn = kn < 1 ? 1 : 1 / kn;
      let kp = Math.pow(C.rotMax / this.rot1p, this.beta);
      kp = kp < 1 ? 1 : 1 / kp;
      const Eup = this.E1p, Eun = this.E1n;
      if (T.load === 1) {
        T.load = 2;
        if (C.stress >= 0) {
          T.rotPu = C.strain - C.stress / (Eup * kp);
          const energy = C.energyD - 0.5 * C.stress / (Eup * kp) * C.stress;
          let damfac = 0;
          if (C.rotMax > this.rot1p) {
            damfac = this.damfc2 * energy / this.energyA;
            damfac += this.damfc1 * (C.rotMax / this.rot1p - 1);
          }
          T.rotMin = C.rotMin * (1 + damfac);
        }
      }
      T.load = 2;
      T.rotMin = T.rotMin < this.rot1n ? T.rotMin : this.rot1n;
      const minmom = this.negEnvStress(T.rotMin);
      const rotlim = this.posEnvRotlim(C.rotMax);
      const rotrel = rotlim < T.rotPu ? rotlim : T.rotPu;
      const rotmp2 = T.rotMin - (1 - this.pinchY) * minmom / (Eun * kn);
      const rotch = rotrel + (rotmp2 - rotrel) * this.pinchX;
      const s = T.strain;
      if (s > T.rotPu) {
        T.tangent = Eup * kp;
        T.stress = C.stress + T.tangent * dStrain;
        if (T.stress <= 0) { T.stress = 0; T.tangent = Eup * 1e-9; }
      } else if (s <= T.rotPu && s > rotch) {
        if (s >= rotrel) { T.stress = 0; T.tangent = Eun * 1e-9; }
        else {
          T.tangent = minmom * this.pinchY / (rotch - rotrel);
          const t1 = C.stress + Eun * kn * dStrain;
          const t2 = (s - rotrel) * T.tangent;
          if (t1 > t2) { T.stress = t1; T.tangent = Eun * kn; } else T.stress = t2;
        }
      } else {
        T.tangent = (1 - this.pinchY) * minmom / (T.rotMin - rotch);
        const t1 = C.stress + Eun * kn * dStrain;
        const t2 = this.pinchY * minmom + (s - rotch) * T.tangent;
        if (t1 > t2) { T.stress = t1; T.tangent = Eun * kn; } else T.stress = t2;
      }
    }
    commit() { cpHysteretic(this.C, this.T); }
    revert() {
      cpHysteretic(this.T, this.C);
      this.eps = this.C.strain; this.sig = this.C.stress; this.tan = this.C.tangent;
    }
  }

  /** Envelope helper: points [[s, e], ...] with positive magnitudes. */
  function envelope(pts) {
    const s = pts.map(p => Math.abs(p[0])), e = pts.map(p => Math.abs(p[1]));
    const E = [];
    let area = 0;
    for (let k = 0; k < s.length; k++) {
      const s0 = k ? s[k - 1] : 0, e0 = k ? e[k - 1] : 0;
      if (!(e[k] > e0)) throw new Error('Invalid args: envelope deformations must increase');
      E.push((s[k] - s0) / (e[k] - e0));
      area += 0.5 * (e[k] - e0) * (s[k] + s0);
    }
    return { s, e, E, n: s.length, area };
  }

  function envStress(P, x) {
    if (x <= 0) return 0;
    for (let k = 0; k < P.n; k++) {
      if (x <= P.e[k]) return k === 0 ? P.E[0] * x : P.s[k - 1] + P.E[k] * (x - P.e[k - 1]);
    }
    const k = P.n - 1;
    if (P.E[k] > 0) return P.s[k] + P.E[k] * (x - P.e[k]);
    return P.s[k];
  }

  function envTangent(P, x) {
    if (x < 0) return P.E[0] * 1e-9;
    for (let k = 0; k < P.n; k++) if (x <= P.e[k]) return P.E[k];
    const k = P.n - 1;
    return P.E[k] > 0 ? P.E[k] : P.E[0] * 1e-9;
  }

  /** Deformation where the softening segment holding x reaches zero force. */
  function envRotlim(P, x) {
    if (x <= P.e[0]) return null;
    let k = P.n - 1;
    for (let i = 1; i < P.n; i++) if (x <= P.e[i]) { k = i; break; }
    if (P.E[k] < 0) return P.e[k - 1] - P.s[k - 1] / P.E[k];
    return null;
  }

  /** Kent-Scott-Park concrete without tension (Concrete01). */
  function cpConcrete01(d, s) { d.minS = s.minS; d.endS = s.endS; d.slope = s.slope; d.strain = s.strain; d.stress = s.stress; d.tangent = s.tangent; return d; }

  class Concrete01 {
    constructor(fpc, epsc0, fpcu, epscu) {
      this.fpc = -Math.abs(fpc); this.epsc0 = -Math.abs(epsc0);
      this.fpcu = -Math.abs(fpcu); this.epscu = -Math.abs(epscu);
      this.Ec0 = 2 * this.fpc / this.epsc0;
      this.E0 = this.Ec0;
      this.C = { minS: 0, endS: 0, slope: this.Ec0, strain: 0, stress: 0, tangent: this.Ec0 };
      this.T = cpConcrete01({}, this.C);
      this.sig = 0; this.tan = this.Ec0; this.eps = 0;
    }
    setTrial(strain) {
      const C = this.C, T = this.T;
      cpConcrete01(T, C);
      const dS = strain - C.strain;
      if (Math.abs(dS) >= EPS) {
        T.strain = strain;
        if (strain > 0) { T.stress = 0; T.tangent = 0; }
        else {
          const temp = C.stress + T.slope * strain - T.slope * C.strain;
          if (strain < C.strain) {
            this.reload();
            if (temp > T.stress) { T.stress = temp; T.tangent = T.slope; }
          } else if (temp <= 0) { T.stress = temp; T.tangent = T.slope; }
          else { T.stress = 0; T.tangent = 0; }
        }
      }
      this.eps = T.strain; this.sig = T.stress; this.tan = T.tangent;
    }
    reload() {
      const T = this.T;
      if (T.strain <= T.minS) {
        T.minS = T.strain;
        this.envelope();
        this.unload();
      } else if (T.strain <= T.endS) {
        T.tangent = T.slope;
        T.stress = T.tangent * (T.strain - T.endS);
      } else { T.stress = 0; T.tangent = 0; }
    }
    envelope() {
      const T = this.T, { fpc, epsc0, fpcu, epscu } = this;
      if (T.strain > epsc0) {
        const eta = T.strain / epsc0;
        T.stress = fpc * (2 * eta - eta * eta);
        T.tangent = this.Ec0 * (1 - eta);
      } else if (T.strain > epscu) {
        T.tangent = (fpc - fpcu) / (epsc0 - epscu);
        T.stress = fpc + T.tangent * (T.strain - epsc0);
      } else { T.stress = fpcu; T.tangent = 0; }
    }
    unload() {
      const T = this.T;
      let tempStrain = T.minS;
      if (tempStrain < this.epscu) tempStrain = this.epscu;
      const eta = tempStrain / this.epsc0;
      let ratio = 0.707 * (eta - 2) + 0.834;
      if (eta < 2) ratio = 0.145 * eta * eta + 0.13 * eta;
      T.endS = ratio * this.epsc0;
      const temp1 = T.minS - T.endS;
      const temp2 = T.stress / this.Ec0;
      if (temp1 > -EPS) T.slope = this.Ec0;
      else if (temp1 <= temp2) { T.endS = T.minS - temp1; T.slope = T.stress / temp1; }
      else { T.endS = T.minS - temp2; T.slope = this.Ec0; }
    }
    commit() { cpConcrete01(this.C, this.T); }
    revert() {
      cpConcrete01(this.T, this.C);
      this.eps = this.C.strain; this.sig = this.C.stress; this.tan = this.C.tangent;
    }
  }

  /** Kent-Scott-Park concrete with linear tension softening (Concrete02). */
  function cpConcrete02(d, s) { d.ecmin = s.ecmin; d.dept = s.dept; d.eps = s.eps; d.sig = s.sig; d.e = s.e; return d; }

  class Concrete02 {
    constructor(fc, epsc0, fcu, epscu, rat, ft, Ets) {
      this.fc = -Math.abs(fc); this.epsc0 = -Math.abs(epsc0);
      this.fcu = -Math.abs(fcu); this.epscu = -Math.abs(epscu);
      Object.assign(this, { rat, ft, Ets });
      this.Ec0 = 2 * this.fc / this.epsc0;
      this.E0 = this.Ec0;
      this.C = { ecmin: 0, dept: 0, eps: 0, sig: 0, e: this.Ec0 };
      this.T = cpConcrete02({}, this.C);
      this.sig = 0; this.tan = this.Ec0; this.eps = 0;
    }
    tens(epsc) {
      const Ec0 = this.Ec0, ft = this.ft, Ets = this.Ets;
      const eps0 = ft / Ec0, epsu = ft * (1 / Ets + 1 / Ec0);
      if (epsc <= eps0) return [epsc * Ec0, Ec0];
      if (epsc <= epsu) return [ft - Ets * (epsc - eps0), -Ets];
      return [0, 1e-10];
    }
    compr(epsc) {
      const { fc, epsc0, fcu, epscu, Ec0 } = this;
      const r = epsc / epsc0;
      if (epsc >= epsc0) return [fc * r * (2 - r), Ec0 * (1 - r)];
      if (epsc > epscu) return [(fcu - fc) * (epsc - epsc0) / (epscu - epsc0) + fc, (fcu - fc) / (epscu - epsc0)];
      return [fcu, 1e-10];
    }
    setTrial(strain) {
      const C = this.C, T = this.T;
      cpConcrete02(T, C);
      const ec0 = this.Ec0;
      const eps = strain, deps = eps - C.eps;
      T.eps = eps;
      if (Math.abs(deps) < EPS) { this.eps = strain; this.sig = T.sig; this.tan = T.e; return; }
      let sig, e;
      if (eps < T.ecmin) {
        [sig, e] = this.compr(eps);
        T.ecmin = eps;
      } else {
        const epsr = (this.fcu - this.rat * ec0 * this.epscu) / (ec0 * (1 - this.rat));
        const sigmr = ec0 * epsr;
        const sigmm = this.compr(T.ecmin)[0];
        const er = (sigmm - sigmr) / (T.ecmin - epsr);
        const ept = T.ecmin - sigmm / er;
        if (eps <= ept) {
          const sigmin = sigmm + er * (eps - T.ecmin);
          const sigmax = er * 0.5 * (eps - ept);
          sig = C.sig + ec0 * deps;
          e = ec0;
          if (sig <= sigmin) { sig = sigmin; e = er; }
          if (sig >= sigmax) { sig = sigmax; e = 0.5 * er; }
        } else {
          const epn = ept + T.dept;
          if (eps <= epn) {
            const sicn = this.tens(T.dept)[0];
            e = T.dept !== 0 ? sicn / T.dept : ec0;
            sig = e * (eps - ept);
          } else {
            const epstmp = eps - ept;
            [sig, e] = this.tens(epstmp);
            T.dept = eps - ept;
          }
        }
      }
      T.sig = sig; T.e = e;
      this.eps = strain; this.sig = sig; this.tan = e;
    }
    commit() { cpConcrete02(this.C, this.T); }
    revert() {
      cpConcrete02(this.T, this.C);
      this.eps = this.C.eps; this.sig = this.C.sig; this.tan = this.C.e;
    }
  }

  /**
   * Concrete with a compression envelope comp(ε ≤ 0) and an optional tension
   * envelope tens(e ≥ 0). Unloading and reloading in compression follow one
   * straight line through the plastic strain; in tension they follow the
   * secant to the crack origin. Used for Concrete04 and Concrete07.
   *   opts.plastic  'karsan' (Karsan-Jirsa, as Concrete01/04) or
   *                 'mander' (Chang-Mander secant rule, as Concrete07)
   *   opts.tension  'zero' (tension strain measured from 0, as Concrete04) or
   *                 'plastic' (measured from the plastic strain, as Concrete07)
   */
  function cpEnvConcrete(d, s) { d.emin = s.emin; d.ep = s.ep; d.etmax = s.etmax; d.dead = s.dead; d.eps = s.eps; d.sig = s.sig; d.tan = s.tan; return d; }

  class EnvConcrete {
    constructor(comp, tens, Ec, epsc0, fc, opts) {
      Object.assign(this, { comp, tens, Ec, epsc0, fc });
      this.opts = Object.assign({ plastic: 'karsan', tension: 'plastic' }, opts || {});
      this.E0 = Ec;
      this.C = { emin: 0, ep: 0, etmax: 0, dead: false, eps: 0, sig: 0, tan: Ec };
      this.T = cpEnvConcrete({}, this.C);
      this.sig = 0; this.tan = Ec; this.eps = 0;
    }
    plastic(emin, smin) {
      let ep;
      if (this.opts.plastic === 'mander') {
        const Esec = this.Ec * (Math.abs(smin) / (this.Ec * Math.abs(this.epsc0)) + 0.57) /
                     (Math.abs(emin / this.epsc0) + 0.57);
        ep = emin - smin / Esec;
      } else {
        const eta = emin / this.epsc0;
        const ratio = eta < 2 ? 0.145 * eta * eta + 0.13 * eta : 0.707 * (eta - 2) + 0.834;
        ep = ratio * this.epsc0;
      }
      if (ep < emin) ep = emin;
      if (ep > 0) ep = 0;
      if (emin - ep < 0 && smin / (emin - ep) > this.Ec) ep = emin - smin / this.Ec;
      return ep;
    }
    setTrial(eps) {
      const C = this.C, T = this.T;
      cpEnvConcrete(T, C);
      T.eps = eps;
      let sig = 0, tan = 0;
      if (T.dead) { sig = 0; tan = 0; }
      else if (eps < T.ep) {
        if (eps <= T.emin) {
          [sig, tan] = this.comp(eps);
          T.emin = eps;
          if (Math.abs(sig) < 1e-9 * Math.abs(this.fc)) T.dead = true;
          else T.ep = this.plastic(eps, sig);
        } else {
          const sm = this.comp(T.emin)[0];
          const d = T.emin - T.ep;
          tan = d < 0 ? sm / d : this.Ec;
          sig = tan * (eps - T.ep);
        }
      } else if (this.tens && !(this.opts.tension === 'zero' && eps <= 0)) {
        const e = this.opts.tension === 'zero' ? eps : eps - T.ep;
        if (e >= T.etmax) {
          [sig, tan] = this.tens(e);
          T.etmax = e;
        } else {
          const st = this.tens(T.etmax)[0];
          tan = T.etmax > 0 ? st / T.etmax : this.Ec;
          sig = tan * e;
        }
      }
      T.sig = sig; T.tan = tan;
      this.eps = eps; this.sig = sig; this.tan = tan;
    }
    commit() { cpEnvConcrete(this.C, this.T); }
    revert() {
      cpEnvConcrete(this.T, this.C);
      this.eps = this.C.eps; this.sig = this.C.sig; this.tan = this.C.tan;
    }
  }

  function concrete04(fc, ec, ecu, Ec, ft, et, beta) {
    fc = -Math.abs(fc); ec = -Math.abs(ec); ecu = -Math.abs(ecu);
    const Esec = fc / ec;
    if (!(Ec > Esec * 1.0001)) throw new Error('Invalid args: Concrete04 needs Ec > fc/epsc0');
    const n = Ec / (Ec - Esec);
    const comp = eps => {
      if (eps < ecu) return [0, 0];
      const x = eps / ec, xn = Math.pow(x, n), D = n - 1 + xn;
      return [fc * x * n / D, Esec * n * (n - 1) * (1 - xn) / (D * D)];
    };
    let tens = null;
    if (ft > 0) {
      const et0 = ft / Ec;
      const b = beta > 0 && beta < 1 ? beta : 0.1;
      tens = e => {
        if (e <= Math.min(et0, et)) return [Ec * e, Ec];
        if (et > et0 && e <= et) {
          const s = ft * Math.pow(b, (e - et0) / (et - et0));
          return [s, s * Math.log(b) / (et - et0)];
        }
        return [0, 0];
      };
    }
    return () => new EnvConcrete(comp, tens, Ec, ec, fc, { plastic: 'karsan', tension: 'zero' });
  }

  /** Tsai curve y(x), dy/dx with a straight descending branch after xcr. */
  function tsai(n, r, xcr) {
    const y = x => n * x / (1 + (n - r / (r - 1)) * x + Math.pow(x, r) / (r - 1));
    const z = x => {
      const D = 1 + (n - r / (r - 1)) * x + Math.pow(x, r) / (r - 1);
      return n * (1 - Math.pow(x, r)) / (D * D);
    };
    const ycr = y(xcr), zcr = z(xcr);
    const xsp = zcr < 0 ? xcr - ycr / zcr : Infinity;
    return x => {
      if (x <= xcr) return [y(x), z(x)];
      if (x <= xsp) return [ycr + zcr * (x - xcr), zcr];
      return [0, 0];
    };
  }

  function concrete07(fc, ec, Ec, ft, et, xp, xn, r) {
    fc = -Math.abs(fc); ec = -Math.abs(ec);
    if (Math.abs(r - 1) < 1e-6 || !(Ec > 0)) throw new Error('Invalid args: Concrete07 needs r != 1 and Ec > 0');
    const cn = tsai(Ec * ec / fc, r, xn);
    const comp = eps => {
      const [y, z] = cn(eps / ec);
      return [fc * y, fc / ec * z];
    };
    let tens = null;
    if (ft > 0 && et > 0) {
      const tn = tsai(Ec * et / ft, r, xp);
      tens = e => {
        const [y, z] = tn(e / et);
        return [ft * y, ft / et * z];
      };
    }
    return () => new EnvConcrete(comp, tens, Ec, ec, fc, { plastic: 'mander', tension: 'plastic' });
  }

  /**
   * Masing-rule hysteresis on a symmetric skeleton f(x ≥ 0) -> [stress, tangent],
   * with memory of open loops. Used for ReinforcingSteel.
   */
  class Masing {
    constructor(f, E0) {
      this.f = f;
      this.E0 = E0;
      this.C = { stack: [], dir: 0, eps: 0, sig: 0, tan: E0 };
      this.T = { stack: [], dir: 0, eps: 0, sig: 0, tan: E0 };
      this.sig = 0; this.tan = E0; this.eps = 0;
    }
    setTrial(eps) {
      const C = this.C, T = this.T;
      const deps = eps - C.eps;
      T.eps = eps;
      if (Math.abs(deps) < EPS) {
        T.stack = C.stack.slice(); T.dir = C.dir; T.sig = C.sig; T.tan = C.tan;
      } else {
        const d = deps > 0 ? 1 : -1;
        const st = C.stack.slice();
        if (C.dir !== 0 && d !== C.dir) st.push([C.eps, C.sig]);
        for (;;) {
          const k = st.length;
          if (!k) break;
          const tgt = k >= 2 ? st[k - 2][0] : -st[0][0];
          if (d * (eps - tgt) >= 0) { st.pop(); if (k >= 2) st.pop(); } else break;
        }
        let s, t;
        if (!st.length) {
          [s, t] = this.f(Math.abs(eps));
          if (eps < 0) s = -s;
        } else {
          const r = st[st.length - 1];
          const [fs, ft] = this.f(d * (eps - r[0]) / 2);
          s = r[1] + d * 2 * fs;
          t = ft;
        }
        T.stack = st; T.dir = d; T.sig = s; T.tan = t;
      }
      this.eps = eps; this.sig = T.sig; this.tan = T.tan;
    }
    commit() {
      const C = this.C, T = this.T;
      C.stack = T.stack.slice(); C.dir = T.dir; C.eps = T.eps; C.sig = T.sig; C.tan = T.tan;
    }
    revert() {
      const C = this.C;
      this.T = { stack: C.stack.slice(), dir: C.dir, eps: C.eps, sig: C.sig, tan: C.tan };
      this.eps = C.eps; this.sig = C.sig; this.tan = C.tan;
    }
  }

  function reinforcingSteel(fy, fu, Es, Esh, esh, eult) {
    const ey = fy / Es;
    if (!(fu > fy && eult > esh && esh >= ey && Esh > 0))
      throw new Error('Invalid args: ReinforcingSteel needs fu > fy, eult > esh >= fy/Es and Esh > 0');
    const p = Esh * (eult - esh) / (fu - fy);
    const tiny = Es * 1e-6, R = 10;
    const f = x => {
      if (x <= esh) {          // Menegotto-Pinto elbow (R = 10) into the yield plateau
        const xi = x / ey, d = 1 + Math.pow(xi, R);
        return [fy * xi / Math.pow(d, 1 / R), Math.max(tiny, Es * Math.pow(d, -1 / R - 1))];
      }
      if (x <= eult) {
        const r = (eult - x) / (eult - esh);
        return [fu + (fy - fu) * Math.pow(r, p), Math.max(tiny, p * (fu - fy) / (eult - esh) * Math.pow(r, p - 1))];
      }
      return [fu, tiny];
    };
    return () => new Masing(f, Es);
  }

  // ------------------------------------------------------------------------
  // Factory: OpenSees material name + argument list (model units) ->
  // function creating a new material instance.
  // ------------------------------------------------------------------------
  const IMPLEMENTED = ['Elastic', 'ENT', 'ElasticPP', 'Steel01', 'Steel02', 'SteelMPF',
    'ReinforcingSteel', 'Hysteretic', 'HystereticSM', 'Concrete01', 'Concrete02',
    'Concrete04', 'Concrete07'];

  function uniaxial(name, args) {
    const bad = need => {
      throw new Error(`WARNING insufficient args / Invalid #args: uniaxialMaterial ${name} ` +
        `needs ${need}, got ${args.length} (${args.join(' ')})`);
    };
    const num = args.every(a => typeof a === 'number' && isFinite(a));
    if (name !== 'HystereticSM' && !num)
      throw new Error(`WARNING invalid type: uniaxialMaterial ${name} arguments must be numbers (${args.join(' ')})`);
    const a = args;
    switch (name) {
      case 'Elastic':
        if (a.length < 1) bad('E');
        return () => new Elastic(a[0]);
      case 'ENT':
        if (a.length < 1) bad('E');
        return () => new ENT(a[0]);
      case 'ElasticPP':
        if (a.length < 2) bad('E epsyP <epsyN eps0>');
        return () => new ElasticPP(a[0], a[1], a.length > 2 ? a[2] : undefined, a[3] || 0);
      case 'Steel01':
        if (a.length !== 3 && a.length !== 7) bad('Fy E0 b <a1 a2 a3 a4>');
        return () => new Steel01(...a);
      case 'Steel02': {
        if (a.length < 3) bad('Fy E0 b <R0 cR1 cR2 <a1 a2 a3 a4 sigInit>>');
        const [Fy, E0, b, R0 = 15, cR1 = 0.925, cR2 = 0.15, a1 = 0, a2 = 1, a3 = 0, a4 = 1, sig0 = 0] = a;
        return () => new Steel02(Fy, Fy, E0, b, b, R0, cR1, cR2, a1, a2, a3, a4, sig0);
      }
      case 'SteelMPF': {
        if (a.length < 8) bad('fyp fyn E0 bp bn R0 cR1 cR2 <a1 a2 a3 a4>');
        const [fyp, fyn, E0, bp, bn, R0, cR1, cR2, a1 = 0, a2 = 1, a3 = 0, a4 = 1] = a;
        return () => new Steel02(fyp, Math.abs(fyn), E0, bp, bn, R0, cR1, cR2, a1, a2, a3, a4);
      }
      case 'ReinforcingSteel':
        if (a.length < 6) bad('fy fu Es Esh esh eult');
        return reinforcingSteel(a[0], a[1], a[2], a[3], a[4], a[5]);
      case 'Hysteretic': {
        let pos, neg, rest;
        if (a.length === 17 || a.length === 18) {
          pos = [[a[0], a[1]], [a[2], a[3]], [a[4], a[5]]];
          neg = [[a[6], a[7]], [a[8], a[9]], [a[10], a[11]]];
          rest = a.slice(12);
        } else if (a.length === 13 || a.length === 14) {
          pos = [[a[0], a[1]], [a[2], a[3]]];
          neg = [[a[4], a[5]], [a[6], a[7]]];
          rest = a.slice(8);
        } else bad('s1p e1p s2p e2p <s3p e3p> s1n e1n s2n e2n <s3n e3n> pinchX pinchY damage1 damage2 <beta>');
        envelope(pos); envelope(neg);
        return () => new Hysteretic(pos, neg, rest[0], rest[1], rest[2], rest[3], rest[4] || 0);
      }
      case 'HystereticSM': {
        const o = { posEnv: [], negEnv: [], pinch: [1, 1], damage: [0, 0], beta: [0] };
        let key = null;
        for (const t of a) {
          if (typeof t === 'string') {
            const k = t.replace(/^-/, '');
            if (!(k in o)) throw new Error(`WARNING invalid type: HystereticSM option ${t} is not supported`);
            key = k;
            o[key] = [];
          } else if (key) o[key].push(t);
          else bad('-posEnv ... -negEnv ...');
        }
        const pairs = v => { const r = []; for (let i = 0; i + 1 < v.length; i += 2) r.push([v[i], v[i + 1]]); return r; };
        const pos = pairs(o.posEnv), neg = o.negEnv.length ? pairs(o.negEnv) : pos;
        for (const [side, p] of [['positive', pos], ['negative', neg]])
          if (p.length < 2 || p.length > 7) throw new Error(`HystereticSM ${side} envelope needs 2-7 points, got ${p.length}`);
        envelope(pos); envelope(neg);
        return () => new Hysteretic(pos, neg, o.pinch[0], o.pinch[1], o.damage[0], o.damage[1], o.beta[0] || 0);
      }
      case 'Concrete01':
        if (a.length < 4) bad('fpc epsc0 fpcu epscu');
        return () => new Concrete01(a[0], a[1], a[2], a[3]);
      case 'Concrete02':
        if (a.length < 7) bad('fpc epsc0 fpcu epscu lambda ft Ets');
        return () => new Concrete02(...a.slice(0, 7));
      case 'Concrete04':
        if (a.length < 4) bad('fc ec ecu Ec <fct et <beta>>');
        return concrete04(a[0], a[1], a[2], a[3], a[4] || 0, a[5] || 0, a.length > 6 ? a[6] : 0.1);
      case 'Concrete07':
        if (a.length < 8) bad('fc ec Ec ft et xp xn r');
        return concrete07(...a.slice(0, 8));
      default:
        throw new Error(`WARNING unknown material type "${name}" - this version supports: ${IMPLEMENTED.join(', ')}`);
    }
  }

  // ========================================================================
  // Part 2 - catalogue (port of materials.py)
  // ========================================================================
  const S = 'steel', C = 'concrete', J = 'spring';

  const P = (name, label, kind, def, extra) => Object.assign({ name, label, kind, default: def, help: '', role: '' }, extra || {});
  const defaultFor = (prm, family) => (family in prm.default ? prm.default[family] : Object.values(prm.default)[0]);

  const fy = d => P('fy', 'Yield strength', 'f', d || { [S]: 420.0, [J]: 150.0 });
  const E_ = (d, name = 'E0') => P(name, 'Initial stiffness', 'k', d || { [S]: 200000.0, [J]: 50000.0 }, { role: 'E' });
  const b_ = d => P('b', 'Strain-hardening ratio b', 'r', d || { [S]: 0.01, [J]: 0.02 });
  const R_ = () => [P('R0', 'R0', 'r', { [S]: 18.0 }), P('cR1', 'cR1', 'r', { [S]: 0.925 }), P('cR2', 'cR2', 'r', { [S]: 0.15 })];
  const HYST_TAIL = () => [
    P('pinchX', 'Pinching x', 'r', { [S]: 1.0, [J]: 0.8 }),
    P('pinchY', 'Pinching y', 'r', { [S]: 1.0, [J]: 0.2 }),
    P('damage1', 'Damage (ductility)', 'r', { [S]: 0.0 }),
    P('damage2', 'Damage (energy)', 'r', { [S]: 0.0 }),
    P('beta', 'Unloading β', 'r', { [S]: 0.0 }, { help: 'Degraded unloading stiffness exponent' })];

  function concBasic(fcu = true, ecu = 0.005) {
    const ps = [P('fc', "Compressive strength f'c", 'f', { [C]: 30.0 }),
                P('epsc0', "Strain at f'c", 'e', { [C]: 0.002 })];
    if (fcu) ps.push(P('fcu', 'Crushing strength', 'f', { [C]: 6.0 }));
    ps.push(P('epscu', 'Ultimate strain', 'e', { [C]: ecu }));
    return ps;
  }

  function confineKentPark(v, K, ecu) {
    const core = Object.assign({}, v, { fc: K * v.fc, epsc0: K * v.epsc0, epscu: ecu });
    if ('fcu' in v) core.fcu = 0.2 * core.fc;
    return core;
  }

  function envArgs(v, s) {
    const env = v.env;
    const pos = env.pos, neg = env.neg || env.pos;
    for (const [side, pts] of [['positive', pos], ['negative', neg]])
      if (!(pts.length >= 2 && pts.length <= 7))
        throw new Error(`HystereticSM ${side} envelope needs 2-7 points, got ${pts.length}`);
    const out = ['-posEnv'];
    for (const [F, e] of pos) out.push(F * s, e);
    out.push('-negEnv');
    for (const [F, e] of neg) out.push(-Math.abs(F) * s, -Math.abs(e));
    return out.concat(['-pinch', v.pinchX, v.pinchY, '-damage', v.damage1, v.damage2, '-beta', v.beta]);
  }

  function parseCustom(text) {
    return String(text).trim().split(/[,\s]+/).filter(Boolean).map(t => {
      const x = Number(t);
      return isNaN(x) ? t : x;
    });
  }

  const LIST = [
    { name: 'Elastic', families: [S, C, J],
      params: [P('E', 'Elastic modulus', 'k', { [S]: 200000.0, [C]: 25000.0, [J]: 50000.0 }, { role: 'E' })],
      args: (v, s) => [v.E * s], E0: v => v.E },
    { name: 'ElasticPP', families: [S, J],
      params: [E_(null, 'E'), P('epsy', 'Yield strain / rotation', 'e', { [S]: 0.0021, [J]: 0.003 })],
      args: (v, s) => [v.E * s, v.epsy], E0: v => v.E, note: 'Elastic-perfectly plastic.' },
    { name: 'ENT', families: [C],
      params: [P('E', 'Elastic modulus', 'k', { [C]: 25000.0 })],
      args: (v, s) => [v.E * s], E0: v => v.E, note: 'Elastic, no tension.' },
    { name: 'Steel01', families: [S, J], params: [fy(), E_(), b_()],
      args: (v, s) => [v.fy * s, v.E0 * s, v.b], E0: v => v.E0,
      note: 'Bilinear with kinematic hardening.' },
    { name: 'Steel02', families: [S, J], params: [fy(), E_(), b_()].concat(R_()),
      args: (v, s) => [v.fy * s, v.E0 * s, v.b, v.R0, v.cR1, v.cR2], E0: v => v.E0,
      note: 'Giuffré-Menegotto-Pinto.' },
    { name: 'SteelMPF', families: [S, J],
      params: [P('fyp', 'Yield strength (+)', 'f', { [S]: 420.0, [J]: 150.0 }),
               P('fyn', 'Yield strength (−), magnitude', 'f', { [S]: 420.0, [J]: 150.0 }),
               E_(), P('bp', 'Hardening ratio (+)', 'r', { [S]: 0.01, [J]: 0.02 }),
               P('bn', 'Hardening ratio (−)', 'r', { [S]: 0.01, [J]: 0.02 }),
               P('R0', 'R0', 'r', { [S]: 20.0 }), P('cR1', 'cR1', 'r', { [S]: 0.925 }),
               P('cR2', 'cR2', 'r', { [S]: 0.15 })],
      args: (v, s) => [v.fyp * s, v.fyn * s, v.E0 * s, v.bp, v.bn, v.R0, v.cR1, v.cR2],
      E0: v => v.E0, note: 'Menegotto-Pinto with asymmetric yield (Kolozvari et al.).' },
    { name: 'ReinforcingSteel', families: [S],
      params: [fy(), P('fu', 'Ultimate strength', 'f', { [S]: 600.0 }), E_(null, 'Es'),
               P('Esh', 'Hardening modulus', 'k', { [S]: 4000.0 }),
               P('esh', 'Strain at hardening onset', 'e', { [S]: 0.01 }),
               P('eult', 'Strain at peak stress', 'e', { [S]: 0.12 })],
      args: (v, s) => [v.fy * s, v.fu * s, v.Es * s, v.Esh * s, v.esh, v.eult], E0: v => v.Es,
      note: 'Chang & Mander backbone with Masing-type cycles (no buckling or fatigue).' },
    { name: 'Hysteretic', families: [S, J],
      params: [P('s1', 'Envelope point 1: force', 'f', { [S]: 420.0, [J]: 150.0 }),
               P('e1', 'Envelope point 1: deformation', 'e', { [S]: 0.0021, [J]: 0.003 }),
               P('s2', 'Envelope point 2: force', 'f', { [S]: 500.0, [J]: 180.0 }),
               P('e2', 'Envelope point 2: deformation', 'e', { [S]: 0.02, [J]: 0.03 }),
               P('s3', 'Envelope point 3: force', 'f', { [S]: 420.0, [J]: 60.0 }),
               P('e3', 'Envelope point 3: deformation', 'e', { [S]: 0.06, [J]: 0.08 })].concat(HYST_TAIL()),
      args: (v, s) => [v.s1 * s, v.e1, v.s2 * s, v.e2, v.s3 * s, v.e3,
                       -v.s1 * s, -v.e1, -v.s2 * s, -v.e2, -v.s3 * s, -v.e3,
                       v.pinchX, v.pinchY, v.damage1, v.damage2, v.beta],
      E0: v => v.s1 / v.e1,
      note: 'Trilinear, symmetric envelope. Use HystereticSM for an asymmetric one.' },
    { name: 'HystereticSM', families: [S, J],
      params: [P('env', 'Backbone envelope', 'env',
                 { [S]: { pos: [[420.0, 0.0021], [500.0, 0.02], [420.0, 0.06]], neg: null },
                   [J]: { pos: [[150.0, 0.003], [180.0, 0.03], [170.0, 0.06], [40.0, 0.1]], neg: null } })]
        .concat(HYST_TAIL()),
      args: envArgs, E0: v => v.env.pos[0][0] / v.env.pos[0][1],
      note: 'Multi-point (2–7) envelope hysteretic model (S. Mazzoni).' },
    { name: 'Concrete01', families: [C], params: concBasic(),
      args: (v, s) => [-v.fc * s, -v.epsc0, -v.fcu * s, -v.epscu],
      E0: v => 4700 * Math.sqrt(v.fc), confine: confineKentPark,
      note: 'Kent-Scott-Park, no tension.' },
    { name: 'Concrete02', families: [C],
      params: concBasic().concat([
        P('lam', 'Unloading slope ratio λ', 'r', { [C]: 0.1 }),
        P('ft', 'Tensile strength', 'f', { [C]: 3.4 }, { help: '≈ 0.62·√f\'c' }),
        P('Ets', 'Tension softening stiffness', 'k', { [C]: 1700.0 }, { help: '≈ ft / 0.002' })]),
      args: (v, s) => [-v.fc * s, -v.epsc0, -v.fcu * s, -v.epscu, v.lam, v.ft * s, v.Ets * s],
      E0: v => 4700 * Math.sqrt(v.fc), confine: confineKentPark,
      note: 'Kent-Scott-Park with linear tension softening.' },
    { name: 'Concrete04', families: [C],
      params: concBasic(false, 0.004).concat([
        P('Ec', 'Initial modulus', 'k', { [C]: 25700.0 }),
        P('ft', 'Tensile strength', 'f', { [C]: 3.4 }),
        P('et', 'Ultimate tensile strain', 'e', { [C]: 0.00013 }),
        P('beta', 'Tension decay β', 'r', { [C]: 0.1 })]),
      args: (v, s) => [-v.fc * s, -v.epsc0, -v.epscu, v.Ec * s, v.ft * s, v.et, v.beta],
      E0: v => v.Ec, confine: confineKentPark, note: 'Popovics.' },
    { name: 'Concrete07', families: [C],
      params: [P('fc', "Compressive strength f'c", 'f', { [C]: 30.0 }),
               P('epsc0', "Strain at f'c", 'e', { [C]: 0.002 }),
               P('Ec', 'Initial modulus', 'k', { [C]: 25700.0 }),
               P('ft', 'Tensile strength', 'f', { [C]: 3.4 }),
               P('et', 'Strain at ft', 'e', { [C]: 0.00013 }),
               P('xp', 'Non-dim. critical tensile strain x⁺', 'r', { [C]: 2.0 }),
               P('xn', 'Non-dim. critical compressive strain x⁻', 'r', { [C]: 2.3 }),
               P('r', 'Shape factor r', 'r', { [C]: 3.87 }, { help: '≈ f\'c / 5.2 − 1.9 (f\'c in MPa)' })],
      args: (v, s) => [-v.fc * s, -v.epsc0, v.Ec * s, v.ft * s, v.et, v.xp, v.xn, v.r],
      E0: v => v.Ec,
      confine: (v, K) => Object.assign({}, v, { fc: K * v.fc, epsc0: K * v.epsc0 }),
      note: 'Chang & Mander envelope (Tsai equation); simplified cyclic rules.' },
    { name: 'Custom', families: [S, C, J],
      params: [P('mat', 'Material name', 'text', { [S]: 'Steel01', [C]: 'Concrete01', [J]: 'Steel01' }),
               P('args', 'Arguments (after the tag)', 'text',
                 { [S]: '420000 200000000 0.01', [C]: '-30000 -0.002 -6000 -0.005', [J]: '150 50000 0.02' },
                 { help: 'Passed exactly as typed, in model units kN and m: stresses in kPa ' +
                         '(MPa × 1000), spring moments in kN·m. Words such as -posEnv are kept as text.' }),
               P('E0', 'Initial stiffness (for elastic properties)', 'k',
                 { [S]: 200000.0, [C]: 25700.0, [J]: 50000.0 },
                 { help: 'Used for the elasticBeamColumn option, the fixity factor and tables.' })],
      args: v => parseCustom(v.args), E0: v => v.E0,
      note: 'Any supported material in OpenSees argument order. The core is not auto-confined.' },
  ];
  const CATALOG = {};
  for (const m of LIST) { m.opsName = m.opsName || m.name; CATALOG[m.name] = m; }

  const DEFAULT = { [S]: 'Steel02', [C]: 'Concrete02', [J]: 'Elastic' };

  const namesFor = family => LIST.filter(m => m.families.includes(family)).map(m => m.name);

  /** [OpenSees material name, argument list] for a material spec. */
  function opsArgs(spec, scale) {
    const m = CATALOG[spec.type];
    const name = spec.type === 'Custom' ? String(spec.v.mat).trim() : m.opsName;
    return [name, m.args(spec.v, scale)];
  }

  /** Function creating new instances of the material described by spec. */
  function maker(spec, scale) {
    const [name, args] = opsArgs(spec, scale);
    return uniaxial(name, args);
  }

  const initialE = spec => CATALOG[spec.type].E0(spec.v);

  function confined(spec, K, epscuCore) {
    const m = CATALOG[spec.type];
    if (!m.confine) return spec;
    return { type: spec.type, v: m.confine(spec.v, K, epscuCore) };
  }

  const isConfinable = spec => !!CATALOG[spec.type].confine;

  function defaultSpec(name, family) {
    const v = {};
    for (const p of CATALOG[name].params) v[p.name] = JSON.parse(JSON.stringify(defaultFor(p, family)));
    return { type: name, v };
  }

  FA.materials = {
    S, C, J, CATALOG, LIST, DEFAULT, IMPLEMENTED, namesFor, opsArgs, maker, uniaxial, initialE,
    confined, isConfinable, defaultSpec, defaultFor, parseCustom,
    models: { Elastic, ENT, ElasticPP, Steel01, Steel02, Hysteretic, Concrete01, Concrete02, EnvConcrete, Masing },
  };
})(globalThis.FA = globalThis.FA || {});
