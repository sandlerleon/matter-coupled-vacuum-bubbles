# -*- coding: utf-8 -*-
"""Nonperturbative test of the granular decay rate: the O(4) bounce allowed to deform around a localized dark-matter source.

The worst case for the decay is a source at the centre of the bubble. A static source is a line along Euclidean time, so with the source on the bubble
axis the field depends on (tau, r) only (r the spatial distance from the source). The bounce is found as the maximum of a constrained minimum:
    S_min(c) = min { S[u] : C[u] = c },     C[u] = int d^4x s(u),     b = max_c S_min(c) - S[u_bg],
with S[u] = int d^4x [ (grad u)^2/2 + p0(u) + theta_s(r) (u - 1)^2/2 ], u_bg the static field around the source (the false vacuum with the source),
and theta_s a uniform sphere of radius r_s carrying the dimensionless coupling q of one particle (the cell-average coupling is q n). C is the
4-volume on the true-vacuum side of the barrier, s(u) = [1 - tanh((u - x_top)/w)]/2 (a linear functional of u would be satisfied cheaply by a small
uniform shift of the false vacuum, so it cannot serve as the collective coordinate). The constraint
is imposed by a stiff quadratic penalty and the minima are found with L-BFGS, continuing in c. For q = 0 the method must reproduce the shooting value
of the zero-density bounce. For q > 0 the increase of the action, db_np(q), is compared with the first-order estimate q I_s used in granular_field.py.

Second part, self-screening: the static field around one source of radius r_s is solved by Newton minimization of the finite-volume energy on a
geometric radial grid, at two resolutions (a case is accepted only if both converge and agree to 1%); the far field is a Yukawa tail whose
amplitude defines an effective coupling q_eff, which saturates at q_eff ~ 4 pi r_s when the source pins the field (as for thin-shell screening of
chameleon fields). Writes nonpert_bounce.json.

    python nonpert_bounce.py
"""
import json
import os
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize
from scipy.integrate import solve_bvp

import vacuum_phase_map as VP

HERE = os.path.dirname(os.path.abspath(__file__))
WJ = json.load(open(os.path.join(HERE, "vacuum_wall.json"), encoding="utf-8"))
ETA = WJ["base"]["eta"]
B0_SHOOT = WJ["walls"]["0.30"]["rows"][0]["b"]
st = sorted(x for x in VP.stationary_tb(ETA, 0.0) if VP.d2p_tb(x, ETA, 0.0) > 0)
X_T0, X_H0 = st[0], st[-1]
M_F = np.sqrt(VP.d2p_tb(X_H0, ETA, 0.0))
X_TOP = sorted(VP.stationary_tb(ETA, 0.0))[1]
p0 = lambda u: (u * u - 1.0) ** 2 / 8.0 + ETA * u / 2.0
dp0 = lambda u: u * (u * u - 1.0) / 2.0 + ETA / 2.0


class Grid:
    def __init__(self, h=0.25, L=22.0, R=22.0, q=0.0, rs=0.5):
        self.h = h
        self.nt, self.nr = int(round(2 * L / h)), int(round(R / h))
        self.tau = -L + (np.arange(self.nt) + 0.5) * h
        self.r = (np.arange(self.nr) + 0.5) * h
        rf = np.arange(self.nr + 1) * h                         # faces
        self.V3 = 4.0 * np.pi / 3.0 * (rf[1:] ** 3 - rf[:-1] ** 3)
        self.A = 4.0 * np.pi * rf[1:] ** 2                      # outer face of cell j
        Vs = 4.0 / 3.0 * np.pi * rs ** 3
        frac = np.clip((np.minimum(rf[1:], rs) ** 3 - np.minimum(rf[:-1], rs) ** 3) / np.maximum(rf[1:] ** 3 - rf[:-1] ** 3, 1e-300), 0.0, 1.0)
        self.th = q / Vs * frac                                 # cell-averaged source density, total q
        self.q, self.rs = q, rs
        self.ubg = self.background()

    def static_energy(self, u):
        du = np.diff(np.append(u, X_H0))
        e = np.sum(self.A * 0.5 * du ** 2 / self.h) + np.sum(self.V3 * (p0(u) + 0.5 * self.th * (u - 1.0) ** 2))
        g = self.V3 * (dp0(u) + self.th * (u - 1.0))
        g[:-1] += -self.A[:-1] * du[:-1] / self.h
        g[1:] += self.A[:-1] * du[:-1] / self.h
        g[-1] += -self.A[-1] * du[-1] / self.h
        return e, g

    def background(self):
        u0 = np.full(self.nr, X_H0)
        res = minimize(self.static_energy, u0, jac=True, method="L-BFGS-B", options=dict(maxiter=20000, gtol=1e-12, ftol=1e-15))
        return res.x

    def action(self, U):
        """S[U] - S[u_bg] and its gradient; U has shape (nt, nr); Dirichlet u_bg at the tau ends and at r = R."""
        h, V3, A = self.h, self.V3, self.A
        ub = self.ubg
        P = np.vstack([ub, U, ub])                               # tau ghosts
        dt = np.diff(P, axis=0)                                  # (nt+1, nr)
        Q = np.hstack([U, np.full((self.nt, 1), X_H0)])          # r ghost
        dr = np.diff(Q, axis=1)                                  # (nt, nr)
        drb = np.diff(np.append(ub, X_H0))
        S = (np.sum(V3 * 0.5 * dt ** 2) / h + np.sum(A * 0.5 * dr ** 2) - self.nt * np.sum(A * 0.5 * drb ** 2)
             + h * np.sum(V3 * (p0(U) + 0.5 * self.th * (U - 1.0) ** 2 - p0(ub) - 0.5 * self.th * (ub - 1.0) ** 2)))
        G = h * V3 * (dp0(U) + self.th * (U - 1.0))
        G += V3 * (dt[:-1] - dt[1:]) / h
        G[:, :-1] += -A[:-1] * dr[:, :-1]
        G[:, 1:] += A[:-1] * dr[:, :-1]
        G[:, -1] += -A[-1] * dr[:, -1]
        return S, G

    def C(self, U, xt=X_TOP, w=0.15):
        th = np.tanh((U - xt) / w)
        return self.h * np.sum(self.V3 * 0.5 * (1.0 - th)), -self.h * self.V3 * 0.5 * (1.0 - th * th) / w


def bounce(q=0.0, rs=0.5, h=0.4, L=20.0, R=20.0, kappa=1e-3, verbose=False):
    t0 = time.time()
    g = Grid(h, L, R, q, rs)
    # starting configuration: the zero-density shooting bounce, scaled to a small size
    T, Rr = np.meshgrid(g.tau, g.r, indexing="ij")
    rho = np.sqrt(T ** 2 + Rr ** 2)
    prof = WJ["profiles"]["0.00"]
    xb = lambda s: np.interp(rho / s, prof["r"], prof["x"], right=X_H0)
    U = np.minimum(g.ubg[None, :], xb(1.0))                     # the zero-density bounce with the source's background where that is larger
    c_b = g.C(g.ubg[None, :] + (xb(1.0) - X_H0))[0]             # constraint value of the zero-density bounce (scale of the scan)
    out = []

    def run_c(c, U):
        def F(x):
            Uu = x.reshape(g.nt, g.nr)
            S, GS = g.action(Uu)
            Cv, GC = g.C(Uu)
            return S + kappa * (Cv - c) ** 2, (GS + 2 * kappa * (Cv - c) * GC).ravel()
        res = minimize(F, U.ravel(), jac=True, method="L-BFGS-B", options=dict(maxiter=20000, maxcor=30, gtol=1e-7, ftol=1e-15))
        U = res.x.reshape(g.nt, g.nr)
        S = g.action(U)[0]
        out.append((float(c / c_b), float(S), float(g.C(U)[0] / c_b), int(res.nit)))
        if verbose:
            print("   c/c_b %.3f  S %.3f  C/c_b %.4f  it %d" % out[-1], flush=True)
        return U

    Us = {}
    Uc = U
    for f in np.arange(1.0, 1.81, 0.1):                          # upward and downward from the zero-density bounce, warm-started
        Uc = run_c(f * c_b, Uc); Us[round(f, 3)] = Uc
    Uc = U
    for f in np.arange(0.9, 0.39, -0.1):
        Uc = run_c(f * c_b, Uc); Us[round(f, 3)] = Uc
    kk = int(np.argmax([o[1] for o in out]))
    f0 = out[kk][0]
    start = Us[min(Us, key=lambda k: abs(k - f0))]
    for f in np.linspace(f0 - 0.08, f0 + 0.08, 9):
        run_c(f * c_b, start)
    out.sort()
    Ss = np.array([o[1] for o in out]); cc = np.array([o[2] for o in out])
    k = int(np.argmax(Ss))
    lo, hi = max(k - 2, 0), min(k + 3, len(Ss))
    pf = np.polyfit(cc[lo:hi], Ss[lo:hi], 2)
    b = float(np.polyval(pf, -pf[1] / (2 * pf[0]))) if pf[0] < 0 else float(Ss[k])
    return dict(q=q, rs=rs, h=h, b=b, b_grid_max=float(Ss[k]), c_peak=float(cc[k]), scan=out, seconds=time.time() - t0,
                ubg_center=float(g.ubg[0]))


def first_order(q, rs, h=0.4, L=20.0, R=20.0):
    """q I_s: the source term of the zero-density bounce (on the same grid) relative to the zero-density false vacuum."""
    g = Grid(h, L, R, q, rs)
    T, Rr = np.meshgrid(g.tau, g.r, indexing="ij")
    prof = WJ["profiles"]["0.00"]
    xb = np.interp(np.sqrt(T ** 2 + Rr ** 2), prof["r"], prof["x"], right=X_H0)
    return float(h * np.sum(g.V3 * 0.5 * g.th * ((xb - 1.0) ** 2 - (X_H0 - 1.0) ** 2)))


def _static_newton(rs, q, n, R=40.0):
    """Static field around one source (uniform sphere of radius rs, coupling q) by Newton minimization of the finite-volume energy on a geometric
    radial grid (tridiagonal Hessian); Dirichlet u = x_h at r = R. Returns cell centres, field, convergence flag."""
    from scipy.linalg import solve_banded
    f = np.concatenate([[0.0], np.geomspace(rs / 50.0, R, n)])
    c = 0.5 * (f[1:] + f[:-1])
    V = 4.0 * np.pi / 3.0 * (f[1:] ** 3 - f[:-1] ** 3)
    Vs = 4.0 / 3.0 * np.pi * rs ** 3
    frac = np.clip((np.minimum(f[1:], rs) ** 3 - np.minimum(f[:-1], rs) ** 3) / (f[1:] ** 3 - f[:-1] ** 3), 0.0, 1.0)
    th = q / Vs * frac
    A = 4.0 * np.pi * f[1:-1] ** 2
    d = np.diff(c)
    AR, dR = 4.0 * np.pi * R ** 2, R - c[-1]
    w = A / d
    E = lambda u: (np.sum(0.5 * w * np.diff(u) ** 2) + 0.5 * AR / dR * (X_H0 - u[-1]) ** 2 + np.sum(V * (p0(u) + 0.5 * th * (u - 1.0) ** 2)))
    u = np.where(c < rs, 1.0, X_H0)
    ok = False
    for it in range(200):
        g = V * (dp0(u) + th * (u - 1.0))
        du = np.diff(u)
        g[:-1] -= w * du
        g[1:] += w * du
        g[-1] -= AR / dR * (X_H0 - u[-1])
        hd = V * ((3.0 * u * u - 1.0) / 2.0 + th)
        hd[:-1] += w
        hd[1:] += w
        hd[-1] += AR / dR
        ab = np.zeros((3, u.size))
        ab[0, 1:] = -w
        ab[1] = hd
        ab[2, :-1] = -w
        step = solve_banded((1, 1), ab, -g)
        e0, t = E(u), 1.0
        while E(u + t * step) > e0 + 1e-4 * t * np.dot(g, step) and t > 1e-8:
            t *= 0.5
        u = u + t * step
        if np.max(np.abs(t * step)) < 1e-13:
            ok = bool(np.all(hd > 0))
            break
    return c, u, ok


def screening(rs, qs, R=40.0, n=2400):
    """Static coupling of one source of radius rs: effective coupling from the Yukawa tail of the nonlinear solution, at two radial resolutions.
    success = Newton converged at both resolutions and q_eff agrees between them to 1%."""
    out = []
    for q in qs:
        res = []
        for nn in (n, 2 * n):
            c, u, ok = _static_newton(rs, q, nn, R)
            rt = np.linspace(4.0 / M_F, 10.0 / M_F, 40)
            Amp = np.median((np.interp(rt, c, u) - X_H0) * rt * np.exp(M_F * rt))
            res.append((float(-4 * np.pi * Amp / (X_H0 - 1.0)), ok, float(u[0])))
        conv = abs(res[1][0] / res[0][0] - 1.0)
        out.append(dict(q=q, q_eff=res[1][0], q_eff_coarse=res[0][0], resolution_change=float(conv), u_center=res[1][2],
                        success=bool(res[0][1] and res[1][1] and conv < 0.01)))
    return out


def _job(a):
    return bounce(*a)


def _cjob(a):
    q, rs, h, L = a
    return bounce(q, rs, h, L, L)


def convergence():
    """Grid and box convergence of the bounce at large coupling: q = 0, 10 and 1e4 (r_s = 0.5) at three grid spacings, and with a larger box."""
    jobs = [(q, 0.5, h, 20.0) for h in (0.4, 0.3, 0.25) for q in (0.0, 10.0, 1e4)] + [(q, 0.5, 0.4, 26.0) for q in (0.0, 10.0, 1e4)]
    with Pool(6) as pool:
        res = pool.map(_cjob, jobs)
    rows = []
    for h, L in ((0.4, 20.0), (0.3, 20.0), (0.25, 20.0), (0.4, 26.0)):
        sel = {r["q"]: r for r, j in zip(res, jobs) if j[2] == h and j[3] == L}
        rows.append(dict(h=h, L=L, b0=sel[0.0]["b"], db_q10=sel[10.0]["b"] - sel[0.0]["b"], db_q1e4=sel[1e4]["b"] - sel[0.0]["b"]))
    out = json.load(open(os.path.join(HERE, "nonpert_bounce.json"), encoding="utf-8"))
    out["convergence"] = rows
    json.dump(out, open(os.path.join(HERE, "nonpert_bounce.json"), "w"), indent=1)
    for r in rows:
        print(r)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "convergence":
        convergence()
        sys.exit()
    if len(sys.argv) > 1 and sys.argv[1] == "saturation":         # saturation at fine resolution: q = 1e3 and 1e4 at h = 0.25 and 0.2
        jobs = [(0.0, 0.5, 0.2, 20.0), (1e3, 0.5, 0.2, 20.0), (1e4, 0.5, 0.2, 20.0), (1e3, 0.5, 0.25, 20.0)]
        with Pool(4) as pool:
            res = pool.map(_cjob, jobs)
        fp = os.path.join(HERE, "nonpert_bounce.json")
        out = json.load(open(fp, encoding="utf-8"))
        b0_25 = [c for c in out["convergence"] if c["h"] == 0.25][0]["b0"]
        out["saturation_fine"] = dict(h02=dict(b0=res[0]["b"], db_q1e3=res[1]["b"] - res[0]["b"], db_q1e4=res[2]["b"] - res[0]["b"]),
                                      h025=dict(b0=b0_25, db_q1e3=res[3]["b"] - b0_25,
                                                db_q1e4=[c for c in out["convergence"] if c["h"] == 0.25][0]["db_q1e4"]))
        json.dump(out, open(fp, "w"), indent=1)
        print(out["saturation_fine"])
        sys.exit()
    if len(sys.argv) > 1 and sys.argv[1] == "screening":          # refresh only the screening block
        fp = os.path.join(HERE, "nonpert_bounce.json")
        out = json.load(open(fp, encoding="utf-8"))
        out["screening"] = {"%.3g" % rs: screening(rs, [0.01, 0.1, 1.0, 10.0, 100.0, 1e3, 1e4]) for rs in (0.01, 0.1, 0.5)}
        json.dump(out, open(fp, "w"), indent=1)
        print({k: all(x["success"] for x in v) for k, v in out["screening"].items()})
        sys.exit()
    t0 = time.time()
    jobs = [(0.0, 0.5), (1.0, 0.5), (3.0, 0.5), (10.0, 0.5), (30.0, 0.5), (100.0, 0.5), (1e3, 0.5), (1e4, 0.5), (1e4, 0.25), (1e4, 1.0)]
    with Pool(5) as pool:
        res = pool.map(_job, jobs)
    b0 = [r for r in res if r["q"] == 0.0][0]["b"]
    for r in res:
        r["db_np"] = r["b"] - b0
        r["db_first_order"] = first_order(r["q"], r["rs"])
        r["ratio_np_over_first"] = r["db_np"] / r["db_first_order"] if r["db_first_order"] else None
    scr = {"%.3g" % rs: screening(rs, [0.01, 0.1, 1.0, 10.0, 100.0, 1e3, 1e4]) for rs in (0.01, 0.1, 0.5)}
    out = dict(eta=ETA, b0_shooting=B0_SHOOT, b0_grid=b0, grid_error=b0 / B0_SHOOT - 1.0, bounces=res, screening=scr, m_f=float(M_F))
    fp = os.path.join(HERE, "nonpert_bounce.json")
    if os.path.exists(fp):                                         # keep the convergence block written by 'python nonpert_bounce.py convergence'
        old = json.load(open(fp, encoding="utf-8"))
        if "convergence" in old:
            out["convergence"] = old["convergence"]
    json.dump(out, open(os.path.join(HERE, "nonpert_bounce.json"), "w"), indent=1)
    print("b0 shooting %.2f, grid %.2f (%.2f%%)" % (B0_SHOOT, b0, 100 * (b0 / B0_SHOOT - 1)))
    for r in res:
        print("q %8g rs %.2f  b %.2f  db_np %8.2f  db_first %10.2f  ratio %s  ubg(0) %.4f  [%.0fs]" % (
            r["q"], r["rs"], r["b"], r["db_np"], r["db_first_order"], ("%.3f" % r["ratio_np_over_first"]) if r["ratio_np_over_first"] else "--",
            r["ubg_center"], r["seconds"]))
    for k, v in scr.items():
        print("screening rs", k, [(x["q"], round(x["q_eff"], 4), round(x["u_center"], 4)) for x in v])
    print("[%.0fs]" % (time.time() - t0))
