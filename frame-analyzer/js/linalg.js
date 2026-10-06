/* Linear algebra for the frame analyzer.
 *
 *   Profile   symmetric skyline (profile) storage with an in-place LDLᵀ
 *             factorisation, solve and matrix-vector product
 *   rcm       reverse Cuthill-McKee ordering (reduces the profile)
 *   symEig    eigenvalues/vectors of a dense symmetric matrix (Householder
 *             tridiagonalisation + implicit QL, after JAMA / EISPACK)
 *   small dense helpers: inverse, Cholesky, least squares
 */
(function (FA) {
  'use strict';

  class SingularError extends Error {
    constructor(eq, msg) {
      super(msg || `singular stiffness matrix: zero pivot at equation ${eq}`);
      this.eq = eq;
    }
  }

  // ------------------------------------------------------------------------
  // Skyline profile. Column j stores rows first[j]..j (upper triangle).
  // The value arrays are kept outside so several matrices (tangent, mass,
  // committed tangent, effective stiffness) can share one profile.
  // ------------------------------------------------------------------------
  class Profile {
    constructor(first) {
      const n = first.length;
      this.n = n;
      this.first = Int32Array.from(first);
      this.diag = new Int32Array(n);
      let p = 0;
      for (let j = 0; j < n; j++) {
        p += j - first[j];
        this.diag[j] = p;
        p += 1;
      }
      this.size = p;
    }

    alloc() { return new Float64Array(this.size); }

    /** a(i,j) += v for i, j in either order. */
    add(a, i, j, v) {
      if (i > j) { const t = i; i = j; j = t; }
      a[this.diag[j] - (j - i)] += v;
    }

    get(a, i, j) {
      if (i > j) { const t = i; i = j; j = t; }
      if (i < this.first[j]) return 0;
      return a[this.diag[j] - (j - i)];
    }

    /** In-place LDLᵀ (stored as unit upper U and D on the diagonal). */
    factor(a) {
      const n = this.n, first = this.first, diag = this.diag;
      for (let j = 0; j < n; j++) {
        const fj = first[j], dj = diag[j];
        for (let i = fj + 1; i < j; i++) {
          const fi = first[i], di = diag[i];
          const m = fi > fj ? fi : fj;
          let s = 0;
          let pi = di - (i - m), pj = dj - (j - m);
          for (let k = m; k < i; k++) s += a[pi++] * a[pj++];
          a[dj - (j - i)] -= s;
        }
        const orig = a[dj];
        let d = orig;
        let pj = dj - (j - fj);
        for (let i = fj; i < j; i++, pj++) {
          const g = a[pj];
          const u = g / a[diag[i]];
          d -= g * u;
          a[pj] = u;
        }
        if (!(Math.abs(d) > 1e-12 * Math.abs(orig)) || !isFinite(d)) throw new SingularError(j);
        a[dj] = d;
      }
    }

    /** Solve with a factored array; b is overwritten with the solution. */
    solve(a, b) {
      const n = this.n, first = this.first, diag = this.diag;
      for (let j = 0; j < n; j++) {
        const fj = first[j];
        let s = b[j];
        let p = diag[j] - (j - fj);
        for (let i = fj; i < j; i++) s -= a[p++] * b[i];
        b[j] = s;
      }
      for (let j = 0; j < n; j++) b[j] /= a[diag[j]];
      for (let j = n - 1; j >= 0; j--) {
        const fj = first[j], x = b[j];
        if (x === 0) continue;
        let p = diag[j] - (j - fj);
        for (let i = fj; i < j; i++) b[i] -= a[p++] * x;
      }
      return b;
    }

    /** y = A x (A unfactored); y is overwritten (or allocated). */
    matvec(a, x, y) {
      const n = this.n, first = this.first, diag = this.diag;
      y = y || new Float64Array(n);
      y.fill(0);
      for (let j = 0; j < n; j++) {
        const fj = first[j], xj = x[j];
        let p = diag[j] - (j - fj);
        let s = 0;
        for (let i = fj; i < j; i++, p++) {
          const v = a[p];
          y[i] += v * xj;
          s += v * x[i];
        }
        y[j] += s + a[p] * xj;
      }
      return y;
    }
  }

  // ------------------------------------------------------------------------
  // Reverse Cuthill-McKee. adj: array of neighbour arrays. Returns perm with
  // perm[new] = old.
  // ------------------------------------------------------------------------
  function rcm(adj) {
    const n = adj.length;
    const deg = adj.map(a => a.length);
    const visited = new Uint8Array(n);
    const order = [];
    while (order.length < n) {
      // start from an unvisited node of minimum degree, then move to a
      // pseudo-peripheral node by one BFS sweep
      let start = -1;
      for (let i = 0; i < n; i++) if (!visited[i] && (start < 0 || deg[i] < deg[start])) start = i;
      start = farthest(adj, start, visited, deg);
      const queue = [start];
      visited[start] = 1;
      for (let h = 0; h < queue.length; h++) {
        const v = queue[h];
        order.push(v);
        const nb = adj[v].filter(w => !visited[w]).sort((p, q) => deg[p] - deg[q]);
        for (const w of nb) { visited[w] = 1; queue.push(w); }
      }
    }
    return order.reverse();
  }

  function farthest(adj, s, visitedGlobal, deg) {
    const seen = new Map([[s, 0]]);
    const q = [s];
    for (let h = 0; h < q.length; h++) {
      const v = q[h];
      for (const w of adj[v]) {
        if (visitedGlobal[w] || seen.has(w)) continue;
        seen.set(w, seen.get(v) + 1);
        q.push(w);
      }
    }
    const maxLevel = seen.get(q[q.length - 1]);
    let best = -1;
    for (const v of q) if (seen.get(v) === maxLevel && (best < 0 || deg[v] < deg[best])) best = v;
    return best;
  }

  // ------------------------------------------------------------------------
  // Dense helpers (row-major Float64Array, n x n)
  // ------------------------------------------------------------------------
  /** Inverse of a small dense matrix by Gauss-Jordan with partial pivoting. */
  function inv(A, n, out) {
    const m = Float64Array.from(A);
    const r = out || new Float64Array(n * n);
    r.fill(0);
    for (let i = 0; i < n; i++) r[i * n + i] = 1;
    for (let c = 0; c < n; c++) {
      let p = c, best = Math.abs(m[c * n + c]);
      for (let i = c + 1; i < n; i++) {
        const v = Math.abs(m[i * n + c]);
        if (v > best) { best = v; p = i; }
      }
      if (!(best > 0) || !isFinite(best)) return null;
      if (p !== c) {
        for (let k = 0; k < n; k++) {
          let t = m[c * n + k]; m[c * n + k] = m[p * n + k]; m[p * n + k] = t;
          t = r[c * n + k]; r[c * n + k] = r[p * n + k]; r[p * n + k] = t;
        }
      }
      const piv = 1 / m[c * n + c];
      for (let k = 0; k < n; k++) { m[c * n + k] *= piv; r[c * n + k] *= piv; }
      for (let i = 0; i < n; i++) {
        if (i === c) continue;
        const f = m[i * n + c];
        if (f === 0) continue;
        for (let k = 0; k < n; k++) {
          m[i * n + k] -= f * m[c * n + k];
          r[i * n + k] -= f * r[c * n + k];
        }
      }
    }
    return r;
  }

  /** Lower Cholesky factor of a dense SPD matrix (row-major); null if not SPD. */
  function cholesky(A, n) {
    const L = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
      let s = A[j * n + j];
      for (let k = 0; k < j; k++) s -= L[j * n + k] * L[j * n + k];
      if (!(s > 0)) return null;
      const d = Math.sqrt(s);
      L[j * n + j] = d;
      for (let i = j + 1; i < n; i++) {
        let t = A[i * n + j];
        for (let k = 0; k < j; k++) t -= L[i * n + k] * L[j * n + k];
        L[i * n + j] = t / d;
      }
    }
    return L;
  }

  /** Least squares min ||A c - b|| for A given as an array of column vectors. */
  function lstsq(cols, b) {
    const k = cols.length;
    const G = new Float64Array(k * k), r = new Float64Array(k);
    for (let i = 0; i < k; i++) {
      r[i] = dot(cols[i], b);
      for (let j = 0; j <= i; j++) G[i * k + j] = G[j * k + i] = dot(cols[i], cols[j]);
    }
    for (let i = 0; i < k; i++) G[i * k + i] *= 1 + 1e-12;
    const Gi = inv(G, k);
    if (!Gi) return new Float64Array(k);
    const c = new Float64Array(k);
    for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) c[i] += Gi[i * k + j] * r[j];
    return c;
  }

  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }

  function norm(a, n) {
    n = n === undefined ? a.length : n;
    let s = 0;
    for (let i = 0; i < n; i++) s += a[i] * a[i];
    return Math.sqrt(s);
  }

  // ------------------------------------------------------------------------
  // Symmetric eigenproblem (dense). A: array of n row arrays (copied).
  // Returns { values (ascending), vectors: rows V[i][j], column j = vector j }.
  // ------------------------------------------------------------------------
  function symEig(A) {
    const n = A.length;
    const V = A.map(r => Float64Array.from(r));
    const d = new Float64Array(n), e = new Float64Array(n);
    tred2(V, d, e, n);
    tql2(V, d, e, n);
    return { values: d, vectors: V };
  }

  function tred2(V, d, e, n) {
    for (let j = 0; j < n; j++) d[j] = V[n - 1][j];
    for (let i = n - 1; i > 0; i--) {
      let scale = 0, h = 0;
      for (let k = 0; k < i; k++) scale += Math.abs(d[k]);
      if (scale === 0) {
        e[i] = d[i - 1];
        for (let j = 0; j < i; j++) {
          d[j] = V[i - 1][j];
          V[i][j] = 0;
          V[j][i] = 0;
        }
      } else {
        for (let k = 0; k < i; k++) {
          d[k] /= scale;
          h += d[k] * d[k];
        }
        let f = d[i - 1];
        let g = Math.sqrt(h);
        if (f > 0) g = -g;
        e[i] = scale * g;
        h -= f * g;
        d[i - 1] = f - g;
        for (let j = 0; j < i; j++) e[j] = 0;
        for (let j = 0; j < i; j++) {
          f = d[j];
          V[j][i] = f;
          g = e[j] + V[j][j] * f;
          for (let k = j + 1; k <= i - 1; k++) {
            g += V[k][j] * d[k];
            e[k] += V[k][j] * f;
          }
          e[j] = g;
        }
        f = 0;
        for (let j = 0; j < i; j++) {
          e[j] /= h;
          f += e[j] * d[j];
        }
        const hh = f / (h + h);
        for (let j = 0; j < i; j++) e[j] -= hh * d[j];
        for (let j = 0; j < i; j++) {
          f = d[j];
          g = e[j];
          for (let k = j; k <= i - 1; k++) V[k][j] -= (f * e[k] + g * d[k]);
          d[j] = V[i - 1][j];
          V[i][j] = 0;
        }
      }
      d[i] = h;
    }
    for (let i = 0; i < n - 1; i++) {
      V[n - 1][i] = V[i][i];
      V[i][i] = 1;
      const h = d[i + 1];
      if (h !== 0) {
        for (let k = 0; k <= i; k++) d[k] = V[k][i + 1] / h;
        for (let j = 0; j <= i; j++) {
          let g = 0;
          for (let k = 0; k <= i; k++) g += V[k][i + 1] * V[k][j];
          for (let k = 0; k <= i; k++) V[k][j] -= g * d[k];
        }
      }
      for (let k = 0; k <= i; k++) V[k][i + 1] = 0;
    }
    for (let j = 0; j < n; j++) {
      d[j] = V[n - 1][j];
      V[n - 1][j] = 0;
    }
    V[n - 1][n - 1] = 1;
    e[0] = 0;
  }

  function tql2(V, d, e, n) {
    for (let i = 1; i < n; i++) e[i - 1] = e[i];
    e[n - 1] = 0;
    let f = 0, tst1 = 0;
    const eps = Math.pow(2, -52);
    for (let l = 0; l < n; l++) {
      tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
      let m = l;
      while (m < n) {
        if (Math.abs(e[m]) <= eps * tst1) break;
        m++;
      }
      if (m > l) {
        let iter = 0;
        do {
          iter++;
          let g = d[l];
          let p = (d[l + 1] - g) / (2 * e[l]);
          let r = Math.hypot(p, 1);
          if (p < 0) r = -r;
          d[l] = e[l] / (p + r);
          d[l + 1] = e[l] * (p + r);
          const dl1 = d[l + 1];
          let h = g - d[l];
          for (let i = l + 2; i < n; i++) d[i] -= h;
          f += h;
          p = d[m];
          let c = 1, c2 = c, c3 = c;
          const el1 = e[l + 1];
          let s = 0, s2 = 0;
          for (let i = m - 1; i >= l; i--) {
            c3 = c2;
            c2 = c;
            s2 = s;
            g = c * e[i];
            h = c * p;
            r = Math.hypot(p, e[i]);
            e[i + 1] = s * r;
            s = e[i] / r;
            c = p / r;
            p = c * d[i] - s * g;
            d[i + 1] = h + s * (c * g + s * d[i]);
            for (let k = 0; k < n; k++) {
              h = V[k][i + 1];
              V[k][i + 1] = s * V[k][i] + c * h;
              V[k][i] = c * V[k][i] - s * h;
            }
          }
          p = -s * s2 * c3 * el1 * e[l] / dl1;
          e[l] = s * p;
          d[l] = c * p;
        } while (Math.abs(e[l]) > eps * tst1 && iter < 300);
      }
      d[l] += f;
      e[l] = 0;
    }
    for (let i = 0; i < n - 1; i++) {
      let k = i, p = d[i];
      for (let j = i + 1; j < n; j++) if (d[j] < p) { k = j; p = d[j]; }
      if (k !== i) {
        d[k] = d[i];
        d[i] = p;
        for (let j = 0; j < n; j++) {
          p = V[j][i]; V[j][i] = V[j][k]; V[j][k] = p;
        }
      }
    }
  }

  /**
   * Generalised symmetric problem K x = λ M x for small dense K, M (row
   * arrays, M SPD). Returns ascending values and M-orthonormal vectors
   * (array of Float64Array, one per mode).
   */
  function genEig(K, M) {
    const n = K.length;
    const Mf = new Float64Array(n * n);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Mf[i * n + j] = M[i][j];
    const L = cholesky(Mf, n);
    if (!L) throw new Error('eigen analysis: mass matrix of the reduced problem is not positive definite');
    const Li = lowerInv(L, n);
    // A = L⁻¹ K L⁻ᵀ
    const T = new Float64Array(n * n);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) {
        let s = 0;
        for (let k = 0; k <= i; k++) s += Li[i * n + k] * K[k][j];
        T[i * n + j] = s;
      }
    const A = [];
    for (let i = 0; i < n; i++) {
      const row = new Float64Array(n);
      for (let j = 0; j < n; j++) {
        let s = 0;
        for (let k = 0; k <= j; k++) s += T[i * n + k] * Li[j * n + k];
        row[j] = s;
      }
      A.push(row);
    }
    for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) A[i][j] = A[j][i] = 0.5 * (A[i][j] + A[j][i]);
    const { values, vectors } = symEig(A);
    const out = [];
    for (let m = 0; m < n; m++) {
      const x = new Float64Array(n);
      // x = L⁻ᵀ y
      for (let i = 0; i < n; i++) {
        let s = 0;
        for (let k = i; k < n; k++) s += Li[k * n + i] * vectors[k][m];
        x[i] = s;
      }
      out.push(x);
    }
    return { values, vectors: out };
  }

  function lowerInv(L, n) {
    const R = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
      R[j * n + j] = 1 / L[j * n + j];
      for (let i = j + 1; i < n; i++) {
        let s = 0;
        for (let k = j; k < i; k++) s += L[i * n + k] * R[k * n + j];
        R[i * n + j] = -s / L[i * n + i];
      }
    }
    return R;
  }

  FA.linalg = { Profile, SingularError, rcm, inv, cholesky, lstsq, dot, norm, symEig, genEig };
})(globalThis.FA = globalThis.FA || {});
