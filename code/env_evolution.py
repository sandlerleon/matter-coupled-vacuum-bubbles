# -*- coding: utf-8 -*-
"""Does the transition order a_tr ~ (1 + delta)^(1/3) survive when the density contrast of a region evolves?

The frozen-contrast relation (viability_extended.py) assumes that a region keeps its density contrast. Voids keep emptying: the linear contrast grows
with the LambdaCDM growth factor D(a), and the nonlinear contrast of a spherical void follows 1 + delta = (1 - delta_L/nu)^(-nu), nu = 1.5. The local
coupling at completion, theta_c = theta_mean(a_tr) (1 + delta(a_tr)), is taken from the frozen calculation (it is the same for all regions up to the
weak dependence of the completion criterion on H), and the completion scale factor is solved again with the evolving contrast. Writes
env_evolution.json.

    python env_evolution.py
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
VE = json.load(open(os.path.join(HERE, "viability_extended.json"), encoding="utf-8"))
OMM = 0.311
NU = 1.5


def growth():
    """Linear growth factor D(a) of LambdaCDM (normalized to D(1) = 1), from the growth ODE in ln a."""
    E2 = lambda a: OMM / a ** 3 + (1 - OMM)

    def rhs(lna, y):
        a = np.exp(lna)
        dlnE = -1.5 * OMM / a ** 3 / E2(a)
        return [y[1], -(2.0 + dlnE) * y[1] + 1.5 * OMM / a ** 3 / E2(a) * y[0]]
    sol = solve_ivp(rhs, (np.log(1e-3), np.log(10.0)), [1e-3, 1e-3], dense_output=True, rtol=1e-10, atol=1e-14)
    D1 = sol.sol(0.0)[0]
    return lambda a: sol.sol(np.log(a))[0] / D1


D = growth()
lin = lambda g: NU * (1.0 - g ** (-1.0 / NU))               # linear contrast for a nonlinear 1 + delta = g
nonlin = lambda dl: (1.0 - dl / NU) ** (-NU)
rows = []
for r in VE["prediction"]:
    if r["mean_f0"] != 3.0 or not r["central"] or r["central"] <= 1.0:
        continue
    g0, a_fr = r["one_plus_delta"], r["central"]
    theta_c = a_fr ** -3 * g0                                  # in units of the mean coupling today
    dl0 = lin(g0)
    f = lambda a: a ** -3 * nonlin(dl0 * D(a)) - theta_c
    a_ev = brentq(f, 1.0, 50.0) if f(1.0) > 0 else 1.0
    rows.append(dict(one_plus_delta_today=g0, a_tr_frozen=a_fr, a_tr_evolving=float(a_ev), one_plus_delta_at_tr=float(nonlin(dl0 * D(a_ev))),
                     shift=float(a_ev / a_fr - 1.0)))
out = dict(nu=NU, D_at=[dict(a=a, D=float(D(a))) for a in (1.0, 1.5, 2.0, 3.0, 5.0)], rows=rows)
json.dump(out, open(os.path.join(HERE, "env_evolution.json"), "w"), indent=1)
print(out["D_at"])
for r in rows:
    print(r)
