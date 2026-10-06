/* Section properties and section models.
 *
 * Section spec (as produced by the interface), dimensions in mm, stresses MPa:
 *   { kind: 'Elastic rectangle', b, h, E }
 *   { kind: 'RC rectangle (fiber)', b, h, cover, d_bar, n_top, n_bot, n_side,
 *     conc: spec, K, epscu_core, crack, steel: spec }
 *   { kind: 'Steel I-section (fiber)', d, bf, tw, tf, steel: spec }
 *
 * Section models work in kN, m. Deformation / force order:
 *   2D: [ε0, κz]            -> [P, Mz]
 *   3D: [ε0, κz, κy, θ]     -> [P, Mz, My, T]
 * A fiber at (y, z) has strain ε = ε0 − y·κz + z·κy, as in OpenSees.
 */
(function (FA) {
  'use strict';
  const M = FA.materials;
  const MPA = 1000.0, MM = 0.001;

  function rectJ(b, h) {
    const a = Math.max(b, h), c = Math.min(b, h);
    return a * c ** 3 * (1 / 3 - 0.21 * c / a * (1 - c ** 4 / (12 * a ** 4)));
  }

  /** Elastic properties in kN, m: E, A, Iz (bending along h), Iy (along b), G, J. */
  function sectionProps(s) {
    const k = s.kind;
    let E, A, Iz, Iy, J, nu;
    if (k === 'Elastic rectangle' || k === 'RC rectangle (fiber)') {
      const b = s.b * MM, h = s.h * MM;
      let f;
      if (k === 'Elastic rectangle') { E = s.E * MPA; f = 1.0; }
      else { E = M.initialE(s.conc) * MPA; f = s.crack; }
      A = b * h; Iz = f * b * h ** 3 / 12; Iy = f * h * b ** 3 / 12; J = f * rectJ(b, h); nu = 0.2;
    } else if (k === 'Steel I-section (fiber)') {
      const d = s.d * MM, bf = s.bf * MM, tw = s.tw * MM, tf = s.tf * MM;
      E = M.initialE(s.steel) * MPA; nu = 0.3;
      A = 2 * bf * tf + (d - 2 * tf) * tw;
      Iz = bf * d ** 3 / 12 - (bf - tw) * (d - 2 * tf) ** 3 / 12;
      Iy = 2 * tf * bf ** 3 / 12 + (d - 2 * tf) * tw ** 3 / 12;
      J = (2 * bf * tf ** 3 + (d - 2 * tf) * tw ** 3) / 3;
    } else throw new Error('Unknown section kind ' + k);
    return { E, A, Iz, Iy, G: E / (2 * (1 + nu)), J };
  }

  function linspace(a, b, n) {
    if (n === 1) return [a];
    const out = [];
    for (let i = 0; i < n; i++) out.push(a + (b - a) * i / (n - 1));
    return out;
  }

  /** Bar centres [y, z] in mm (y along the depth h). Cover is to the bar centre. */
  function rcBarLayout(s) {
    const yb = s.h / 2 - s.cover, zb = s.b / 2 - s.cover;
    const bars = [];
    for (const z of linspace(zb, -zb, s.n_top)) bars.push([yb, z]);
    for (const z of linspace(zb, -zb, s.n_bot)) bars.push([-yb, z]);
    const side = linspace(-yb, yb, s.n_side + 2).slice(1, -1);
    for (const y of side) { bars.push([y, zb]); bars.push([y, -zb]); }
    return bars;
  }

  /** [As_total mm², ρ] of an RC rectangle. */
  function rcSteelSummary(s) {
    const aBar = Math.PI * s.d_bar ** 2 / 4;
    const As = aBar * rcBarLayout(s).length;
    return [As, As / (s.b * s.h)];
  }

  // ------------------------------------------------------------------------
  // Section models
  // ------------------------------------------------------------------------
  class ElasticSection {
    constructor(dim, sp) {
      this.order = dim === 2 ? 2 : 4;
      const n = this.order;
      const diag = dim === 2 ? [sp.E * sp.A, sp.E * sp.Iz] : [sp.E * sp.A, sp.E * sp.Iz, sp.E * sp.Iy, sp.G * sp.J];
      this.ks = new Float64Array(n * n);
      for (let i = 0; i < n; i++) this.ks[i * n + i] = diag[i];
      this.diag = diag;
      this.e = new Float64Array(n);
      this.Ce = new Float64Array(n);
      this.s = new Float64Array(n);
    }
    setTrial(e) {
      for (let i = 0; i < this.order; i++) { this.e[i] = e[i]; this.s[i] = this.diag[i] * e[i]; }
    }
    commit() { this.Ce.set(this.e); }
    revert() { this.setTrial(this.Ce); }
  }

  class FiberSection {
    /** fibers: array of [y, z, A, makeMaterial]; GJ used in 3D. */
    constructor(dim, fibers, GJ) {
      this.dim = dim;
      this.order = dim === 2 ? 2 : 4;
      const nf = fibers.length;
      this.nf = nf;
      this.y = new Float64Array(nf);
      this.z = new Float64Array(nf);
      this.A = new Float64Array(nf);
      this.mats = new Array(nf);
      let a = 0, qy = 0, qz = 0;
      fibers.forEach((f, i) => { a += f[2]; qy += f[0] * f[2]; qz += f[1] * f[2]; });
      const yBar = qy / a, zBar = qz / a;   // geometric centroid, as OpenSees
      fibers.forEach((f, i) => {
        this.y[i] = f[0] - yBar;
        this.z[i] = f[1] - zBar;
        this.A[i] = f[2];
        this.mats[i] = f[3]();
      });
      this.GJ = GJ || 0;
      const n = this.order;
      this.e = new Float64Array(n);
      this.Ce = new Float64Array(n);
      this.s = new Float64Array(n);
      this.ks = new Float64Array(n * n);
      this.sumUp();
    }
    setTrial(e) {
      const { y, z, mats, nf } = this;
      this.e.set(e);
      if (this.dim === 2) {
        const e0 = e[0], kz = e[1];
        for (let i = 0; i < nf; i++) mats[i].setTrial(e0 - y[i] * kz);
      } else {
        const e0 = e[0], kz = e[1], ky = e[2];
        for (let i = 0; i < nf; i++) mats[i].setTrial(e0 - y[i] * kz + z[i] * ky);
      }
      this.sumUp();
    }
    sumUp() {
      const { y, z, A, mats, nf, s, ks } = this;
      if (this.dim === 2) {
        let P = 0, Mz = 0, k00 = 0, k01 = 0, k11 = 0;
        for (let i = 0; i < nf; i++) {
          const m = mats[i], a = A[i], yi = y[i];
          const fa = m.sig * a, ea = m.tan * a;
          P += fa; Mz -= fa * yi;
          k00 += ea; k01 -= ea * yi; k11 += ea * yi * yi;
        }
        s[0] = P; s[1] = Mz;
        ks[0] = k00; ks[1] = k01; ks[2] = k01; ks[3] = k11;
      } else {
        let P = 0, Mz = 0, My = 0, k00 = 0, k01 = 0, k02 = 0, k11 = 0, k12 = 0, k22 = 0;
        for (let i = 0; i < nf; i++) {
          const m = mats[i], a = A[i], yi = y[i], zi = z[i];
          const fa = m.sig * a, ea = m.tan * a;
          P += fa; Mz -= fa * yi; My += fa * zi;
          k00 += ea; k01 -= ea * yi; k02 += ea * zi;
          k11 += ea * yi * yi; k12 -= ea * yi * zi; k22 += ea * zi * zi;
        }
        s[0] = P; s[1] = Mz; s[2] = My; s[3] = this.GJ * this.e[3];
        ks.fill(0);
        ks[0] = k00; ks[1] = k01; ks[2] = k02;
        ks[4] = k01; ks[5] = k11; ks[6] = k12;
        ks[8] = k02; ks[9] = k12; ks[10] = k22;
        ks[15] = this.GJ;
      }
    }
    commit() {
      for (const m of this.mats) m.commit();
      this.Ce.set(this.e);
    }
    revert() {
      for (const m of this.mats) m.revert();
      this.e.set(this.Ce);
      this.sumUp();
    }
  }

  /** Rectangular patch -> fibers (OpenSees "patch rect"). */
  function patchRect(out, mk, nfy, nfz, yI, zI, yJ, zJ) {
    const dy = (yJ - yI) / nfy, dz = (zJ - zI) / nfz;
    for (let i = 0; i < nfy; i++)
      for (let j = 0; j < nfz; j++)
        out.push([yI + (i + 0.5) * dy, zI + (j + 0.5) * dz, Math.abs(dy * dz), mk]);
  }

  /**
   * Returns a function creating a new section model for spec s
   * (each integration point needs its own material state).
   */
  function sectionMaker(s, dim) {
    const sp = sectionProps(s);
    if (s.kind === 'Elastic rectangle') return () => new ElasticSection(dim, sp);
    const fibers = [];
    if (s.kind === 'RC rectangle (fiber)') {
      const cover = M.maker(s.conc, MPA);
      const core = M.maker(M.confined(s.conc, s.K, s.epscu_core), MPA);
      const steel = M.maker(s.steel, MPA);
      const b = s.b * MM, h = s.h * MM, c = s.cover * MM;
      if (!(b > 0 && h > 0 && c > 0 && 2 * c < Math.min(b, h)))
        throw new Error('Invalid section: check b, h and the cover (0 < 2·cover < min(b, h))');
      const y1 = h / 2, z1 = b / 2, nf = 16;
      const nz = dim === 2 ? 1 : nf;
      patchRect(fibers, core, nf, nz, -y1 + c, -z1 + c, y1 - c, z1 - c);
      patchRect(fibers, cover, nf, 1, -y1, -z1, y1, -z1 + c);
      patchRect(fibers, cover, nf, 1, -y1, z1 - c, y1, z1);
      patchRect(fibers, cover, 2, nz, -y1, -z1 + c, -y1 + c, z1 - c);
      patchRect(fibers, cover, 2, nz, y1 - c, -z1 + c, y1, z1 - c);
      const As = Math.PI * (s.d_bar * MM) ** 2 / 4;
      for (const [y, z] of rcBarLayout(s)) fibers.push([y * MM, z * MM, As, steel]);
    } else if (s.kind === 'Steel I-section (fiber)') {
      const st = M.maker(s.steel, MPA);
      const d = s.d * MM, bf = s.bf * MM, tw = s.tw * MM, tf = s.tf * MM;
      if (!(d > 2 * tf && bf >= tw && tw > 0 && tf > 0))
        throw new Error('Invalid section: check d, bf, tw, tf (d > 2·tf, bf ≥ tw > 0)');
      const [nzf, nzw] = dim === 2 ? [1, 1] : [8, 2];
      patchRect(fibers, st, 4, nzf, d / 2 - tf, -bf / 2, d / 2, bf / 2);
      patchRect(fibers, st, 4, nzf, -d / 2, -bf / 2, -d / 2 + tf, bf / 2);
      patchRect(fibers, st, 16, nzw, -d / 2 + tf, -tw / 2, d / 2 - tf, tw / 2);
    } else throw new Error('Unknown section kind ' + s.kind);
    const GJ = dim === 3 ? sp.G * sp.J : 0;
    return () => new FiberSection(dim, fibers, GJ);
  }

  FA.sections = { sectionProps, rcBarLayout, rcSteelSummary, sectionMaker, ElasticSection, FiberSection, MPA, MM };
})(globalThis.FA = globalThis.FA || {});
