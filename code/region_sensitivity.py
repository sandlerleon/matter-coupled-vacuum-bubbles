# -*- coding: utf-8 -*-
"""Sensitivity of the transitioned-region results to the energy given to the dark matter at the wall, and of the isolated-bubble trapping result to
the exterior cosmology.

Region evolution (bubble_gravity.Daughter): the number of dark-matter quanta at the crossing multiplied by 1 (adiabatic) or (1 + r^2)/2r (sudden, the
1+1-dimensional simulation of wave_wall_1d.py gives 97% of it), with the energy taken from the wall. Part of the sudden energy is kinetic (the quanta
leave the wall at (r^2-1)/(r^2+1)), so the two runs bracket the evolution.

Isolated bubble (bubble_handoff.isolated_consistency): the trapping test is repeated with the sudden crossing and with three exterior cosmologies (the
observed one with Omega_Lambda shifted by +-0.02, and an exterior whose dark matter carries the held-vacuum coupling energy). The test identifies a
trapped surface (rho_in R^2/3 > 1, an apparent horizon) within the thin-wall junction model; a global event horizon is not established.
Writes region_sensitivity.json.

    python region_sensitivity.py
"""
import json
import os

import numpy as np

import bubble_gravity as BG
import bubble_handoff as BH
import matter_coupling as MC

HERE = os.path.dirname(os.path.abspath(__file__))
WJ = json.load(open(os.path.join(HERE, "vacuum_wall.json"), encoding="utf-8"))
TJ = json.load(open(os.path.join(HERE, "interface_transport.json"), encoding="utf-8"))
WW = json.load(open(os.path.join(HERE, "wave_wall_1d.json"), encoding="utf-8"))
base, TW = WJ["base"], WJ["trigger_wall"]
th_deg = WJ["walls"]["0.30"]["theta_deg"]
r = TJ["dark_matter"]["A_true"] / TJ["dark_matter"]["A_false"]
SUDDEN = (1 + r * r) / (2 * r)


def ev(res, name):
    e = [x for x in res["events"] if x["event"] == name]
    return e[0] if e else None


out = {"sudden_jump": SUDDEN, "regions": {}, "isolated": {}}
for lab, jump in (("adiabatic", 1.0), ("sudden", SUDDEN)):
    D = BG.Daughter(base["eta"], base["f0"] * th_deg, base["q"], TW["theta"], dm_jump=jump)
    rp = D.run("percolated", n_cycles=1)
    nuc, ta, rs = ev(rp, "nucleation"), ev(rp, "turnaround"), ev(rp, "restoration")
    out["regions"][lab] = dict(a_max_over_a_nuc=ta["a_max_over_a_nuc"] if ta else None, lna_turn=ta["lna"] if ta else None,
                               lna_restoration=rs["lna"] if rs else None, dm_energy_from_wall=nuc.get("dm_energy_from_wall"),
                               dV_eff=nuc.get("dV_eff"))
    D = BG.Daughter(base["eta"], base["f0"] * th_deg, base["q"], TW["theta"], dm_jump=jump)
    ri = D.run("isolated", n_cycles=1)
    exteriors = {"observed": None}
    if lab == "adiabatic":
        exteriors.update({"Omega_Lambda +0.02": 1.02 / 1.0, "Omega_Lambda -0.02": 0.98 / 1.0, "held coupling energy in the exterior": "coupling"})
    for ext, mod in exteriors.items():
        saved = (MC.VC, MC.RHO_DM0)
        if isinstance(mod, float):
            MC.VC = saved[0] * (1.0 + (mod - 1.0) / 0.6889)
        elif mod == "coupling":
            MC.RHO_DM0 = saved[1] * TJ["dark_matter"]["A_false"]
        try:
            ic = BH.isolated_consistency({"isolated": ri})
        finally:
            MC.VC, MC.RHO_DM0 = saved
        out["isolated"]["%s, %s" % (lab, ext)] = {k: ic[k] for k in ("t_nucleation", "t_turnaround", "t_trapped", "t_interior_denser", "t_restoration",
                                                                     "t_bounce")}
json.dump(out, open(os.path.join(HERE, "region_sensitivity.json"), "w"), indent=1)
print("sudden jump %.3f" % SUDDEN)
for k, v in out["regions"].items():
    print(k, v)
for k, v in out["isolated"].items():
    print(k, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()})
