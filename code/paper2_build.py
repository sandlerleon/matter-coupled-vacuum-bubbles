# -*- coding: utf-8 -*-
"""Paper 2 (Physics of the Dark Universe): density-dependent vacuum bubble nucleation in matter-coupled dark-sector cosmology.

    python paper2_build.py     -> manuscript/PDU_Vacuum_Bubbles.{docx,tex}, supplement, cover letter
                                  (after interface_validation.py, granularity.py, make_int_figures.py, make_gran_figure.py)
"""
import datetime
import json
import os
import platform
import re
import sys

import numpy as np
import scipy
from docx import Document
from docx.shared import Pt, Inches

import copyedit as CE
import int_blocks as IB
import render as RN

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(REPO, "manuscript")
FIGDIR = os.path.join(REPO, "figures")
J = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
if not IB.VJ["all_passed"]:
    sys.exit("interface checks have not all passed; run interface_validation.py")
DOIS = J("dois.json")
COMP = J("companion.json") if os.path.exists(os.path.join(HERE, "companion.json")) else {}
GITHUB = "https://github.com/sandlerleon/matter-coupled-vacuum-bubbles"
REFS = dict(J("refs_all.json"))
REFS.update({
    "Companion": "L. Sandler, Matter-dilution-triggered recollapse and conditional nonsingular cycles in a metastable dark-energy model. Preprint, "
                 "Zenodo (2026). https://doi.org/%s" % COMP.get("paper", "(to be assigned)"),
    "HOTW17": "L. Hui, J.P. Ostriker, S. Tremaine, E. Witten, Ultralight scalars as cosmological dark matter. Phys. Rev. D 95, 043541 (2017). "
              "https://doi.org/10.1103/PhysRevD.95.043541",
    "PWW83": "J. Preskill, M.B. Wise, F. Wilczek, Cosmology of the invisible axion. Phys. Lett. B 120, 127--132 (1983). "
             "https://doi.org/10.1016/0370-2693(83)90637-8",
    "TG79": "S. Tremaine, J.E. Gunn, Dynamical role of light neutral leptons in cosmology. Phys. Rev. Lett. 42, 407--410 (1979). "
            "https://doi.org/10.1103/PhysRevLett.42.407",
    "GMW14": "R. Gregory, I.G. Moss, B. Withers, Black holes as bubble nucleation sites. J. High Energy Phys. 03, 081 (2014). "
             "https://doi.org/10.1007/JHEP03(2014)081",
    "BGM15": "P. Burda, R. Gregory, I.G. Moss, Vacuum metastability with black holes. J. High Energy Phys. 08, 114 (2015). "
             "https://doi.org/10.1007/JHEP08(2015)114",
    "SvdW04": "R.K. Sheth, R. van de Weygaert, A hierarchy of voids: much ado about nothing. Mon. Not. R. Astron. Soc. 350, 517--538 (2004). "
              "https://doi.org/10.1111/j.1365-2966.2004.07661.x",
})
SV, DBJ, GRJ = J("survival.json"), J("density_bounce.json"), J("granularity.json")
TITLE = ("Density-Dependent Vacuum Bubble Nucleation in Matter-Coupled Dark-Sector Cosmology: Granularity of the Source, Wall Dynamics, Transmission, "
         "and Gravitational Viability")
texsci, sci = IB.texsci, IB.sci
PM, WJ, TJ, GJ, CJ, EJ, VJ, HJ, VM = IB.PM, IB.WJ, IB.TJ, IB.GJ, IB.CJ, IB.EJ, IB.VJ, IB.HJ, IB.VM
COH, PS, TW, DM, LNG, NOREG, DW = IB.COH, IB.PS, IB.TW, IB.DM, IB.LNG, IB.NOREG, IB.DW
PC, WB, IS, RP = HJ["percolation"], HJ["wall_bounds"], HJ["isolated"], HJ["reproduction"]
VAR = {v["variant"]: v for v in HJ["variants"]}
W0 = VAR["wall energy as massive scalar quanta until the bounce (w = 0)"]
W13 = VAR["wall energy as dark radiation until the bounce (w = 1/3)"]
TM = VM["tunneling"]
B_P, B_I = IB.B_P, IB.B_I
ev = IB.ev
TB = WJ["walls"]["0.30"]["rows"]
r0 = TB[0]
r6 = [r for r in TB if abs(r["theta_over_deg"] - 0.6) < 1e-9][0]
deg = WJ["walls"]["0.30"]["degenerate"]
VC = IB.VC
GT, G0 = GRJ["trigger"], GRJ["today"]
VOID = {round(v["one_plus_delta"], 2): v for v in GRJ["voids"]}
DV_ = 0.1                                       # density of the emptiest regions relative to the mean (deepest voids)
F0_LOC = WJ["base"]["f0"]                        # the baseline theta_0/theta_deg, read as the local value in those regions
F0_MEAN = F0_LOC / DV_
THR_T = 1.0 / (F0_MEAN * WJ["a_trigger"] ** -3)  # 1 + delta above which the held vacuum is still the true vacuum at the trigger
NEARC = GRJ["near_eta_c_max"]
NEAR01 = max((r for r in GRJ["near_eta_c"] if r["theta_over_deg"] >= 0.1), key=lambda r: r["ratio"])
VCOMP = GRJ["void_completion"]
HIGHLIGHTS = ["Dilution makes a matter-held vacuum metastable below a degeneracy density",
              "Particle dark matter cannot trigger it: the voids between particles decay first",
              "Wave dark matter between 3 meV and 0.25 eV restores the smooth mean field",
              "Decay starts in the emptiest voids; galactic haloes keep the old vacuum",
              "Walls run away, become transparent to dark matter and carry the released energy"]

B = []
add = B.append
ABSTRACT = (
    "Can a declining dark-matter density switch on the decay of a vacuum, and what do the resulting vacuum walls do? We study a scalar field whose "
    "coupling to dark matter is centred on a false vacuum, so that the matter term holds that vacuum while the density is high. In the mean-field "
    "description dilution turns the held vacuum from the true into a metastable one at a degeneracy density, exactly $\\theta=\\eta/2$ in the natural "
    "units of the potential, and Coleman--De Luccia tunneling switches on below it. We show that this description has a narrow domain of validity. For "
    "a vacuum at the dark-energy scale the field responds within its Compton length, about %.2f mm, inside which the dark-matter energy is only "
    "%.2f$m_\\varphi^4$. Particle dark matter of any mass is then a set of point sources between which the field sees zero density, and the vacuum would "
    "have decayed long ago ($e^{%.0f}$ expected bubbles in our past light cone). The mean field is restored only for wave dark matter that is heavier "
    "than the field and still classical, %.0f meV $\\lesssim m_{\\rm dm}\\lesssim$ %.2f eV. For such dark matter the transition is environmental: it begins "
    "in the emptiest voids, protecting the past requires a mean coupling today close to the degeneracy density, and overdense regions, including all "
    "galactic haloes, keep the held vacuum. Within the voids the wall is about %.1f mm thick, runs away, becomes transparent to the dark matter (a "
    "fraction $%s$ is reflected) and carries the released energy; the transition completes void by void within a few hundredths of a Hubble time, and "
    "the transitioned volume recollapses in its negative vacuum, whereas an isolated Hubble-scale bubble ends as a black hole. A continuation through a "
    "bounce would need a thermal relic and is not supported by the dark matter that the trigger requires."
    % (GT["compton_length_mm"], GT["rho_over_mphi4"], GRJ["particle_dm"]["ln_expected_bubbles_past_light_cone"], GT["m_phi_meV"],
       GT["classical_max_mass_meV"] / 1e3, TW["thickness_mm"], texsci(DM["filtered_fraction"], 0)))
add(("title", TITLE))
add(("author", "Leon Sandler\nIndependent Researcher, Northbrook, Illinois, USA\nsandler.leon@gmail.com"))
add(("abstract", ABSTRACT))
add(("keywords", "vacuum decay; bubble nucleation; coupled dark matter; wave dark matter; dark energy; cosmic voids; Israel junction"))

# ============================================================================================ 1
add(("h1", "1. Introduction"))
add(("p", "False-vacuum decay proceeds by the nucleation of bubbles of the true vacuum, whose rate follows from the Euclidean bounce \\cite{Coleman,CC} and "
          "whose gravitational corrections were derived by Coleman and De Luccia \\cite{CDL}. The theory is usually applied to a vacuum whose stability does not "
          "depend on its surroundings. In a universe filled with dark matter, however, a scalar that couples to the dark matter feels an effective potential "
          "that depends on the dark-matter density, as in chameleon and symmetron screening \\cite{KW,HK,HKLM11} and in coupled dark energy \\cite{Amendola}. "
          "If the dark energy is the energy of such a scalar held in a vacuum by the matter term, the dilution of the dark matter by the expansion changes "
          "the vacuum structure itself, and the decay of the vacuum is placed under the control of the cosmic density. This paper asks three questions. Can a "
          "declining dark-matter density turn a stable vacuum into a metastable one? When is the density that the field feels the smooth mean density "
          "that such an argument assumes? And what do the resulting vacuum walls do: what crosses them, how do they move, and what becomes of the region they "
          "enclose?"))
add(("p", "Environment-dependent vacua are familiar from screened modified gravity, where the field's minimum and mass track the local density "
          "\\cite{KW,HK,HKLM11}; in those applications the field's Compton length is long compared with the separation of the matter particles, and the "
          "density can be averaged. Vacuum decay out of a metastable state has been studied for vacua whose stability is fixed, from the Standard Model Higgs "
          "\\cite{IRS,EGR} to landscapes of many vacua \\cite{BP,FKMS}, including the signatures of bubble collisions \\cite{Kleban,AJ}. Combining the two "
          "raises a question that neither literature had to face: a vacuum at the dark-energy scale has a Compton length below a millimetre, and whether the "
          "dark matter looks smooth on that scale depends on what the dark matter is. The interfaces add a further question. A wall between vacua in which the "
          "dark-matter mass differs is a mass step; slow walls reflect particles whose momentum is below the step, which has been proposed as a way to filter "
          "dark matter at first-order transitions \\cite{BKL20}, while fast walls transmit them and the pressure they exert is bounded \\cite{BM09,BM17}. "
          "Bubbles of a negative vacuum enclose regions that recollapse, and bubbles in an expanding exterior can be cut off behind horizons "
          "\\cite{BGG,FG,FGG}."))
add(("p", "The contributions are the following. (i) The phase structure of a matter-coupled potential is mapped as a function of density and tilt "
          "(Section 2). (ii) The conditions under which the field feels a smooth density are derived; they exclude particle dark matter of any mass and leave "
          "a window for wave dark matter heavier than the field, a bound that follows from the dark-energy scale itself (Section 3). (iii) The bounce, the "
          "wall and the density-dependent action are computed and validated, and the nucleation history is followed in the emptiest regions, where it begins, "
          "with the consequences of the environment for the protection of the past and for the regions that keep the held vacuum (Section 4). (iv) The "
          "transmission of fields and of the dark matter through the wall is computed (Section 5), (v) the gravitational dynamics of the wall and the "
          "completion of the transition are computed with explicit bounds on the approximations (Section 6), and (vi) the region enclosed by the wall is "
          "evolved as a cosmology (Section 7). (vii) Controls that remove the coupling, a viability map and a parameter survey identify where the conditions "
          "hold (Section 8)."))
add(("p", "Units are $\\hbar=c=k_B=1$, $M_{\\rm Pl}$ is the reduced Planck mass, and every number quoted is read from the output of the script named in "
          "parentheses in the archived code. The paper does not require a cosmological bounce. A companion paper treats the symmetron version of the trigger, "
          "whose release is a classical roll on cosmological scales, and its continuation into cycles \\cite{Companion}; the two papers share the class of "
          "models but not their results."))

# ============================================================================================ 2
add(("h1", "2. Matter-coupled vacuum metastability"))
add(("p", "A scalar coupled conformally to dark matter has the action"))
add(("eq", r"S=\int d^4x\sqrt{-g}\left[\frac{M_{\rm Pl}^2R}{2}-\frac{1}{2}(\nabla\varphi)^2-V_0(\varphi)\right]+S_{\rm dm}\left[A^2(\varphi)g_{\mu\nu},\psi_{\rm dm}\right]"))
add(("p", "and, where the dark matter can be described by a smooth nonrelativistic density, the conserved Einstein-frame density $\\rho_*$ enters the field "
          "equation through $V_{\\rm eff}(\\varphi,\\rho_*)=V_0(\\varphi)+\\rho_*A(\\varphi)$ \\cite{Amendola,HK}. When this mean-field description holds is "
          "the subject of Section 3; here it is assumed. We take the tilted quartic $V_0$ in the dimensionless form"))
add(("eq", r"p(x;\eta,\theta)=\frac{(x^2-1)^2}{8}+\frac{\eta x}{2}+\frac{\theta}{2}(x-1)^2,\qquad x=\frac{\varphi}{f},\qquad \theta=\frac{\rho_*}{\lambda f^2M^2},\qquad V_{\rm eff}=\lambda f^4\,p+{\rm const}"))
add(("p", "where $\\eta$ tilts the potential so that $x\\approx-1$ is the zero-density true vacuum and $x\\approx1$ the false one, and the coupling "
          "$A(\\varphi)=1+(\\varphi-f)^2/2M^2$ is centred on the false vacuum, so that the matter term $\\theta(x-1)^2/2$ vanishes there and deepens it. "
          "Because $\\rho_*\\propto a^{-3}$, the mean $\\theta=\\theta_0a^{-3}$ decreases monotonically as the universe expands (vacuum_phase_map.py)."))
add(("p", "The stationary points solve $x^3+(2\\theta-1)x+(\\eta-2\\theta)=0$. Two minima coexist for $\\theta_{\\rm low}<\\theta<\\theta_{\\rm high}$, "
          "with $\\theta_{\\rm low}=0$ for $\\eta<\\eta_c=2/(3\\sqrt3)$, and the window edges are roots of the cubic discriminant (Fig. 1b). The matter-held "
          "minimum $x_h\\approx1$ is the true vacuum above the degeneracy density, which is exactly"))
add(("eq", r"\theta_{\rm deg}=\frac{\eta}{2},\qquad p(x;\eta,\eta/2)=\frac{(x^2-1)^2}{8}+\frac{\eta}{4}\left(x^2+1\right)"))
add(("p", "because $p$ is then even in $x$. Below $\\theta_{\\rm deg}$ the held vacuum is metastable: dilution turns a stable vacuum into a decaying one, "
          "and tunneling is meaningful for $0<\\theta<\\eta/2$. Above $\\theta_{\\rm high}$ (%.3f for $\\eta=0.30$) the other minimum disappears, so at high "
          "density the matter term alone determines the vacuum. For $\\eta>\\eta_c$ the held minimum itself disappears below $\\theta_{\\rm low}$ (for example "
          "%.3f at $\\eta=0.45$), and the release is a classical roll rather than a tunneling event (Fig. 1a)."
          % (PM["tunneling_profiles"]["0.30"]["window"][1], PM["tunneling_profiles"]["0.45"]["window"][0])))
add(("p", "The same structure is absent in the symmetron form of the coupling, $A=1+\\varphi^2/2M^2$ with a symmetric double well: there the held point "
          "$\\varphi=0$ loses stability exactly when the broken minima appear, the occupied minimum is never metastable, and Coleman--De Luccia tunneling "
          "does not apply. Its interfaces are domain walls between $\\pm v$, which form if the roll is not seeded coherently; they would carry up to %.0f%% of "
          "the critical density per Hubble volume, and a tilt $\\epsilon\\gtrsim%s\\,\\mu^2v$ suffices to prevent them. The coupling energy of the dark matter "
          "between the vacua is $\\rho_*\\Delta A$ exactly; the first-order estimate $\\delta(\\Delta V)=\\rho_m\\Delta\\ln m$ agrees with it to %.2f%% when the "
          "coupling changes the rate by a factor $e$ and to %.0f%% when it moves the action across the survival window. The tunneling branch is the one in "
          "which density-dependent bubble nucleation can occur, and the rest of the paper concerns it."
          % (100 * COH["domain_walls"][-1]["wall_energy_per_Hubble_volume_over_rho_crit"], texsci(COH["eps_hat_min"], 0),
             100 * (min(c["exact_over_log"] for c in PM["coupling_check"]) - 1), 100 * (max(c["exact_over_log"] for c in PM["coupling_check"]) - 1))))
add(("fig", "fig1_vacuum_landscape.png", "Fig. 1. (a) Effective potential at five densities ($\\eta=0.30$): dilution turns the matter-held minimum near "
            "$x=1$ from the true into the false vacuum at $\\theta_{\\rm deg}=\\eta/2$. (b) Phase diagram: number of minima and which is lower, with the "
            "degeneracy line and the lines where a minimum disappears; the gray line marks $\\eta_c$.", 6.3))

# ============================================================================================ 3 granularity
add(("h1", "3. When does the field feel a smooth dark-matter density?"))
add(("h2", "3.1 The scales of the problem"))
add(("p", "The effective potential $V_0+\\rho_*A$ is a mean-field statement. The field responds to sources within its Compton length "
          "$m_\\varphi^{-1}$, with $m_\\varphi^2=\\lambda f^2\\,\\partial_x^2p(x_h)$ at the held vacuum, and on the time scale $m_\\varphi^{-1}$; the bounce "
          "itself has a radius of a few Compton lengths. Normalized to the observed dark energy (Section 4.3), $m_\\varphi=%.2f$ meV at the trigger and %.2f meV "
          "today, a Compton length of %.3f mm (granularity.py). The mean dark-matter energy density at the trigger is %.2f meV$^4$, so"
          % (GT["m_phi_meV"], G0["m_phi_meV"], GT["compton_length_mm"], GT["rho_dm_meV4"])))
add(("eq", r"\frac{\rho_{\rm dm}}{m_\varphi^4}=%.3f\ \ (\text{trigger}),\qquad %.3f\ \ (\text{today})" % (GT["rho_over_mphi4"], G0["rho_over_mphi4"])))
add(("p", "Less than one field-mass worth of dark matter lies in a Compton volume. Whether the source is nevertheless smooth depends on what the dark "
          "matter is, and three cases must be distinguished (Table 1, Fig. 2a)."))
add(("h2", "3.2 Particle dark matter"))
sp = GRJ["spacing_over_compton"]
add(("p", "If the dark matter consists of particles whose occupation of phase space is small, the source of the field is a set of points. The mean field "
          "applies only if many of them lie inside a Compton volume, $n_{\\rm dm}\\gg m_\\varphi^3$, which requires $m_{\\rm dm}\\ll\\rho_{\\rm dm}/m_\\varphi^3"
          "=%.2f$ meV. A particle description, on the other hand, needs an occupation per de Broglie cell below one, $m_{\\rm dm}\\gg"
          "(\\rho_{\\rm dm}/v^3)^{1/4}=%.0f$ meV for a velocity dispersion $v=10^{-3}$. The two conditions exclude each other, and fermions lighter than a few "
          "hundred eV are excluded in any case by phase-space bounds \\cite{TG79}. For thermal relics the mismatch is large: the mean separation of 1 MeV "
          "particles is %.0f times the Compton length and that of 1 GeV particles %.0f times. Each particle is a screened point source around which the field "
          "is pinned within a Compton length, as for the thin-shell objects of chameleon screening \\cite{KW}; between them the field sees zero density. The "
          "held vacuum then decays in the space between the particles at the zero-density rate, $\\ln(\\Gamma/H^4)=%.0f$ in the baseline, and the expected "
          "number of bubbles in our past light cone is $e^{%.0f}$. Particle dark matter cannot trigger a density-dependent decay of a vacuum at this scale; with "
          "it the vacuum would have decayed long ago. The bound is least severe where it does not matter: in haloes the density is highest, but there the "
          "held vacuum is the true vacuum and nothing decays (Section 4.4); in the voids, where the decay must begin, the density is lower than the mean and "
          "the bound is stronger."
          % (GT["particle_max_mass_meV"], GT["classical_max_mass_meV"], sp["1 MeV"]["over_compton"], sp["1 GeV"]["over_compton"],
             GRJ["particle_dm"]["ln_gamma_over_H4_zero_density"], GRJ["particle_dm"]["ln_expected_bubbles_past_light_cone"])))
add(("h2", "3.3 Wave dark matter"))
add(("p", "A bosonic dark matter with a large occupation number is a classical field oscillating at its mass \\cite{HOTW17}, and its coupling through "
          "$A^2(\\varphi)g_{\\mu\\nu}$ gives the field a source proportional to the instantaneous trace of its stress-energy. Two regimes follow. If the dark "
          "matter is lighter than the field, $m_{\\rm dm}\\ll m_\\varphi$, the field follows the oscillation; the source passes through zero twice per period, "
          "and tunneling, which takes a time of order $m_\\varphi^{-1}$, proceeds in those intervals at the zero-density rate, with the same consequence as for "
          "particles. If it is heavier, $m_{\\rm dm}\\gg m_\\varphi$, the source oscillates faster than the field can respond, and the field feels its time "
          "average, which for a nonrelativistic field is $\\rho_*A(\\varphi)$. The source is also smooth in space, because the field is coherent over its de "
          "Broglie length $1/(m_{\\rm dm}v)$, which exceeds the Compton length for $m_{\\rm dm}\\ll m_\\varphi/v=%.1f$ eV. Together with classicality, the mean "
          "field of Section 2 applies for"
          % (GT["coherence_max_mass_meV"] / 1e3)))
add(("eq", r"m_\varphi=%.1f\ {\rm meV}\ \ll\ m_{\rm dm}\ \ll\ \left(\frac{\rho_{\rm dm}}{v^3}\right)^{1/4}=%.0f\ {\rm meV}" % (GRJ["wave_window_meV"][0], GRJ["wave_window_meV"][1])))
add(("p", "at the trigger, a window of about two decades, which widens slightly in the emptier, colder voids where the trigger begins. Dark matter in this "
          "mass range is a nonthermal condensate, produced for instance by misalignment \\cite{PWW83}, and behaves as cold dark matter on all scales that "
          "have been observed \\cite{HOTW17}. The trigger of this paper therefore specifies the dark matter: a classical bosonic field of meV to sub-eV mass. "
          "Everything that follows assumes it."))
add(("p", "The averaging can be made explicit. The conformal coupling gives the field equation the source $-(A'/A)\\,T_{\\rm dm}$, with $T_{\\rm dm}$ the "
          "trace of the dark-matter stress-energy. For a homogeneous nonrelativistic field oscillating at $m_{\\rm dm}$ the pressure is "
          "$p=-\\rho\\cos(2m_{\\rm dm}t)$, so $T_{\\rm dm}=-\\rho\\,[1+3\\cos(2m_{\\rm dm}t)]$: the average is the $-\\rho$ of a pressureless fluid, and the "
          "oscillating part drives the field at the frequency $2m_{\\rm dm}$. The response of a field of mass $m_\\varphi$ to that drive is smaller than its "
          "response to a static source by $m_\\varphi^2/(4m_{\\rm dm}^2-m_\\varphi^2)$, so the residual modulation of the density the field feels is "
          "$3(m_\\varphi/2m_{\\rm dm})^2$, below 1%% for $m_{\\rm dm}\\ge10m_\\varphi$ and growing to order one as $m_{\\rm dm}$ approaches $m_\\varphi$, "
          "where the two regimes meet. Fluctuations of the dark-matter density on scales between the de Broglie length and the size of a void are part of "
          "the local density of Section 4.4, not of its granularity."))
rowsT1 = [["particle, any mass", "point sources; $n\\gg m_\\varphi^3$ needs $m<%.2f$ meV, particle regime needs $m>%.0f$ meV" % (GT["particle_max_mass_meV"],
                                                                                                                         GT["classical_max_mass_meV"]),
           "zero-density rate between particles", "excluded"],
          ["fermion, light", "phase-space (Tremaine--Gunn) bound", "--", "excluded"],
          ["wave, $m_{\\rm dm}\\ll m_\\varphi$", "field follows the oscillation", "zero-density rate at the nodes", "excluded"],
          ["wave, $m_\\varphi\\ll m_{\\rm dm}\\ll%.0f$ meV" % GT["classical_max_mass_meV"], "time average, coherent over the de Broglie length",
           "mean field $\\rho_*A(\\varphi)$, local density", "viable"]]
add(("table", "Table 1. Coupled dark matter and the density the field feels at the trigger ($m_\\varphi=%.1f$ meV, $\\rho_{\\rm dm}=%.2f$ meV$^4$, $v=10^{-3}$)."
              % (GT["m_phi_meV"], GT["rho_dm_meV4"]), ["Dark matter", "Condition", "What the field feels", "Trigger"], rowsT1, None))
add(("h2", "3.4 Why the bound is generic"))
add(("p", "The smallness of $\\rho_{\\rm dm}/m_\\varphi^4$ is not a property of the baseline. With $\\lambda f^4=(1+q)V_f/\\Delta\\tilde V_0$ fixed by the "
          "observed dark energy and $m_\\varphi^4=\\lambda^2f^4\\,(\\partial_x^2p)^2$,"))
add(("eq", r"\frac{\rho_{\rm dm}}{m_\varphi^4}=\frac{\rho_{\rm dm}}{V_f}\,\frac{\Delta\tilde V_0}{1+q}\,\frac{1}{\lambda\,(\partial_x^2p)^2},\qquad \lambda=\frac{b(\eta,\theta)}{B}"))
add(("p", "The first factor is of order one because the dark matter and the dark energy have comparable densities, which is what makes a dilution "
          "trigger possible in the first place; and $\\lambda$ is fixed by the requirement that the transition happen within a Hubble time, $B\\approx%.0f$. "
          "Over all %d combinations of tilt, density and $q$ in the bounce tables the ratio stays below %.1f. Toward $\\eta_c$, where the held minimum becomes "
          "shallow, both $\\lambda$ and $\\partial_x^2p$ fall and the ratio grows, but only at vanishing density: it reaches %.0f at $\\eta=%.3f$ and "
          "$\\theta=0$, and stays below %.1f once $\\theta\\ge0.1\\,\\theta_{\\rm deg}$, the densities at which the decay is triggered. Particle dark matter "
          "would need a ratio above $1/v\\approx10^3$. The coincidence of the dark-matter and dark-energy densities, which motivates a dilution trigger, is "
          "thus what forbids it for particle dark matter."
          % (TW["B"], GRJ["structural_scan_n"], GRJ["structural_scan_max"]["ratio"], NEARC["ratio"], NEARC["eta"], NEAR01["ratio"])))
add(("fig", "fig2_granularity.png", "Fig. 2. (a) Conditions for the field to feel a smooth dark-matter density at the trigger: particle dark matter "
            "(red) needs two incompatible conditions; wave dark matter (blue) is averaged where it is heavier than the field and still classical (green "
            "band). (b) Expected number of bubbles in the past light cone of a region whose density is a fraction $1+\\delta$ of the mean, for three mean "
            "couplings today ($\\eta=0.30$, $q=1$, $B_0=%.0f$): the emptiest regions decide whether the past is protected." % WJ["base"]["B0"], 6.6))

# ============================================================================================ 4 nucleation
add(("h1", "4. Density-dependent bubble nucleation"))
add(("h2", "4.1 The bounce and its action"))
add(("p", "Where the mean field applies, the O(4)-symmetric bounce solves, in units of $f$ for the field and $(\\sqrt\\lambda f)^{-1}$ for lengths,"))
add(("eq", r"x''+\frac{3}{r}x'=\partial_xp(x;\eta,\theta),\qquad x'(0)=0,\qquad x(\infty)=x_h,\qquad b(\eta,\theta)=2\pi^2\int_0^\infty r^3\left[\frac{x'^2}{2}+p(x)-p(x_h)\right]dr"))
add(("p", "and the action is $B=b/\\lambda$ \\cite{Coleman,CC}. Gravity enters only through $g=f^2/M_{\\rm Pl}^2$, of order $10^{-60}$ at the dark-energy "
          "scale, so the flat-space bounce is exact for all practical purposes \\cite{CDL,Parke}; Section 6 confirms that the wall's own gravity is negligible. "
          "The bounce is found by overshoot and undershoot on the central value with an eighth-order Runge--Kutta integrator; where the central value lies "
          "within $10^{-5}$ of the true minimum (the thin-wall regime) shooting becomes exponentially sensitive and collocation on the boundary-value problem "
          "is used instead (vacuum_wall.py). At $\\theta=0$ the actions reproduce an independent table to $%s$, and where both methods apply they agree to the "
          "tolerance of the collocation."
          % texsci(max(abs(v["b"] / v["table"] - 1) for v in DBJ["check_theta0"].values()) + 1e-16, 0)))
add(("h2", "4.2 The wall"))
add(("p", "The interface follows from the field equations rather than being imposed. At the degeneracy density the static planar wall, $x''=\\partial_xp$, has "
          "the first integral $z(x)=\\int dx/\\sqrt{2(p-p_{\\rm min})}$, tension $\\sigma=\\int\\sqrt{2(p-p_{\\rm min})}\\,dx=%.3f\\,\\sqrt\\lambda f^3$ and 10--90%% "
          "thickness $%.1f/(\\sqrt\\lambda f)$ for $\\eta=0.30$; the symmetric case $\\eta=\\theta=0$ reproduces $\\sigma=2/3$ and $x=\\tanh(z/2)$. Below the "
          "degeneracy density a static planar wall does not exist, since the energy difference accelerates it, and the wall is read from the bounce: its tension "
          "$\\sigma_b=\\int x'^2dr$ and its thickness change little with density, while the radius grows from $%.1f$ to $%.1f$ and the action from $b=%.0f$ to "
          "$%.0f$ between $\\theta=0$ and $0.6\\,\\theta_{\\rm deg}$ (Table 2, Fig. 3a,b). The thin-wall action"
          % (deg["sigma"], deg["thickness"], r0["R_w"], r6["R_w"], r0["b"], r6["b"])))
add(("eq", r"B_{\rm tw}=\frac{27\pi^2\sigma^4}{2\,\Delta V^3},\qquad \sigma=\int_{x_t}^{x_h}\sqrt{2\,[p(x)-p_{\rm lin}(x)]}\,dx"))
add(("p", "with $p_{\\rm lin}$ the straight line through the two minima, which removes the energy difference from the tension, is recovered within 1%% from "
          "$0.6\\,\\theta_{\\rm deg}$ upwards for $\\eta=0.20$ and 0.30, and overestimates the thick-wall bounce near zero density ($b/B_{\\rm tw}=%.2f$ at "
          "$\\theta=0$ for $\\eta=0.30$). The density dependence of the action is therefore dominated by $\\Delta V(\\theta)$, which falls linearly to zero at "
          "$\\theta_{\\rm deg}$, while the tension is nearly constant." % r0["b_over_B_tw"]))
rows2 = [["%.2f" % r["theta_over_deg"], "%.3f" % r["dV"], "%.0f" % r["b"], "%.2f" % r["b_over_B_tw"], "%.1f" % r["R_w"], "%.2f" % r["thickness"],
          "%.3f" % r["sigma_bounce"], r["method"]] for r in TB if abs(r["theta_over_deg"] * 10 - round(r["theta_over_deg"] * 10)) < 1e-6]
add(("table", "Table 2. Bounce and wall against density ($\\eta=0.30$, units $\\lambda f^4$, $(\\sqrt\\lambda f)^{-1}$ and $\\sqrt\\lambda f^3$).",
     ["$\\theta/\\theta_{\\rm deg}$", "$\\Delta V$", "$b$", "$b/B_{\\rm tw}$", "$R_w$", "thickness", "$\\sigma_b$", "method"], rows2, None))
add(("h2", "4.3 Physical scale, survival and the trigger"))
add(("p", "We normalize to the observed universe: the false vacuum carries the observed dark-energy density $V_f$, the true vacuum $V_t=-qV_f$ with $q=1$, "
          "so that $\\lambda f^4=(1+q)V_f/\\Delta\\tilde V_0$, and the action at zero density is $B_0=b(0)/\\lambda=%.0f$ for $\\eta=0.30$. Then $f=%.2f$ meV and the "
          "wall scale $1/\\sqrt\\lambda f$ is %.3f mm. The rate per unit four-volume is $\\Gamma=\\mathcal A\\,e^{-B(\\theta)}$ with the one-loop prefactor "
          "estimated as $(B/2\\pi)^2$ times the fourth power of the false-vacuum mass scale \\cite{CC}, and the relevant comparison is with the Hubble "
          "four-volume, $\\Gamma/H^4$. The probability that no bubble has nucleated in the past light cone of a present observer is $\\exp(-\\int\\Gamma\\,dV_4)$; "
          "in the observed $\\Lambda$CDM history the four-volume of the past light cone is $V_4=%.4f\\,H_\\Lambda^{-4}$, %.0f times smaller than the Euclidean "
          "four-sphere often used as a shortcut, so that survival requires only $\\Gamma/H_\\Lambda^4\\lesssim%.0f$ for a constant rate \\cite{GW,TWW}. With the "
          "coupling the rate in the past was suppressed by the larger density, and the survival integral is evaluated with the density history "
          "$\\theta=\\theta_0(1+z)^3$. In the future the rate rises as the matter dilutes (Fig. 3c), and the transition completes, in the sense that bubbles "
          "percolate, once $\\Gamma/H^4>9/4\\pi$ \\cite{GW}. The prefactor is the least certain input: a full one-loop determinant \\cite{DM} could "
          "differ from the estimate by a factor comparable to $(B/2\\pi)^2$ itself, $e^{%.1f}$ at the trigger. Because the rate grows at $\\beta=%.0f\\,H$ "
          "(Section 6), such a factor moves the completion by %.2f Hubble times, and it moves the boundary of the region with a protected past by the same "
          "%.1f units in $B_0$; the qualitative results do not depend on it."
          % (WJ["base"]["B0"], PS["f_meV"], PS["length_unit_mm"], SV["V4_HL4"], SV["ratio_euclid_to_lorentzian"], SV["gamma_over_hl4_max_95"],
             2 * np.log(TW["B"] / (2 * np.pi)), PC["beta_over_H"], 2 * np.log(TW["B"] / (2 * np.pi)) / PC["beta_over_H"], 2 * np.log(TW["B"] / (2 * np.pi)))))
add(("p", "Because the field responds to the local density (Section 3.3), these statements hold region by region, with $\\theta_0$ the local value. For "
          "a region with $\\theta_0=%.1f\\,\\theta_{\\rm deg}$, the baseline of this paper, the past light cone is clear and the transition completes at "
          "$a=%.2f$, when $\\theta=%.3f\\,\\theta_{\\rm deg}$ and $B=%.0f$; the nucleated bubble has radius %.2f mm, wall thickness %.2f mm and tension "
          "%.0f meV$^3$."
          % (F0_LOC, WJ["a_trigger"], TW["theta_over_deg"], TW["B"], TW["R_w_mm"], TW["thickness_mm"], TW["sigma_meV3"])))
add(("h2", "4.4 The environment: voids first, haloes last"))
add(("p", "The rate depends exponentially on the local density, so the emptiest regions decide both whether the past is protected and when the decay "
          "begins. The deepest cosmic voids have central densities of 10--20%% of the mean \\cite{SvdW04}. If a region with $1+\\delta$ times the mean "
          "density is followed through the whole past light cone with that contrast (conservative, since voids were less empty in the past), the past light "
          "cone of a region at the mean density is clear for any mean coupling above $\\theta_0=%.2f\\,\\theta_{\\rm deg}$, but a region at %.0f%% or %.0f%% "
          "of the mean requires $\\theta_0\\ge%.2f$ or $%.2f\\,\\theta_{\\rm deg}$ on average (Fig. 2b; granularity.py). With the mean coupling of Section 4.3 "
          "taken as the mean rather than the void value, regions at %.0f%% of the mean or less would already have decayed. A protected past therefore requires a "
          "mean coupling today close to the degeneracy density, and the baseline of Section 4.3 is to be read as the local value in voids at %.0f%% of the "
          "mean, $\\theta_0=%.1f\\,\\theta_{\\rm deg}$ on average."
          % (VOID[1.0]["min_f0_mean"], 100 * 0.2, 100 * 0.1, VOID[0.2]["min_f0_mean"], VOID[0.1]["min_f0_mean"],
             100 * max(k for k, v in VOID.items() if not v["baseline_past_clear"]), 100 * DV_, F0_MEAN)))
add(("p", "With that mean, the held vacuum is still the true vacuum wherever $1+\\delta>\\theta_{\\rm deg}/\\bar\\theta(a)$: at the trigger this is any region "
          "denser than %.2f times the mean, and the threshold rises as $a^3$. A bubble nucleated in a void expands through all connected regions in which the "
          "held vacuum is metastable, and its wall stalls where the local density reaches the degeneracy value, because the driving pressure, proportional to "
          "$\\Delta V_{\\rm eff}$, vanishes there and changes sign beyond. The transition therefore advances from the voids into sheets and filaments as the "
          "universe expands, following the contour $\\theta_{\\rm loc}=\\theta_{\\rm deg}$, while virialized haloes, whose density does not dilute, keep the "
          "held vacuum indefinitely: a galactic halo at about two hundred times the mean density has $\\theta\\approx%.0f\\,\\theta_{\\rm deg}$, far above "
          "$\\theta_{\\rm high}=%.2f\\,\\theta_{\\rm deg}$, where only the held minimum exists. The dark energy measured inside galaxies is therefore unchanged "
          "by the transition, and the haloes become islands of the held vacuum bounded by walls whose tension energy is negligible (%s kg per square metre "
          "for the baseline tension). Separate voids are separated by such walls at the trigger and must each nucleate a bubble; since the rate grows at "
          "$\\beta=%.0f\\,H$ (Section 6), a void of comoving radius 5--50 Mpc does so within %.3f--%.3f Hubble times of the global completion criterion."
          % (THR_T, 200 * F0_MEAN, GRJ["environment"]["theta_high_over_deg"],
             texsci(TW["sigma_meV3"] * 1e-9 / (1.97327e-7) ** 2 * 1.602e-19 / (2.998e8) ** 2, 0).join("$$"),
             PC["beta_over_H"], VCOMP["50 Mpc"]["delay_hubble_times"], VCOMP["5 Mpc"]["delay_hubble_times"])))
fr = lambda a: a ** 3 / F0_MEAN
add(("p", "The advance of the transition can be read off the threshold $1+\\delta_*(a)=a^3\\theta_{\\rm deg}/\\bar\\theta_0$: %.2f at the trigger, %.1f at "
          "$a=2$, %.1f at $a=3$ and %.0f at $a=5$. Sheets and the outskirts of filaments, a few times denser than the mean, are swept within a Hubble time of the "
          "trigger; the cores of filaments, ten or more times the mean, only later; virialized structures never, because their density does not fall. Two "
          "features distinguish this from vacuum decay seeded by inhomogeneities elsewhere in the literature. Compact objects such as black holes catalyse "
          "the decay of a metastable vacuum by lowering the action in their vicinity \\cite{GMW14,BGM15}; here the seed is the absence of matter, and the "
          "decay is suppressed, not catalysed, where matter has collapsed. And the wall does not need to overcome the denser regions: it stops where the "
          "held vacuum becomes the true vacuum, so that the outcome of the transition is a map of the cosmic density field at the time of the trigger."
          % (fr(WJ["a_trigger"]), fr(2.0), fr(3.0), fr(5.0))))
add(("fig", "fig3_interface_properties.png", "Fig. 3. (a) Wall tension (solid) and thickness/10 (dashed) of the O(4) bounce against density for three tilts. "
            "(b) Bounce action $b=\\lambda B$, numerical (points) and thin-wall (dotted). (c) $\\ln(\\Gamma/H^4)$ along the dilution history of a region with "
            "$\\theta_0=%.1f\\,\\theta_{\\rm deg}$, with the completion threshold and the trigger." % F0_LOC, 6.6))

# ============================================================================================ 5
add(("h1", "5. Interface tension and field transmission"))
sp_ = TJ["species"]
g1 = sp_["1.00"]
add(("p", "A wall makes the mass of every field coupled to $\\varphi$ depend on position. For a test scalar with $\\mathcal L_{\\rm int}=-\\frac12g^2\\varphi^2\\psi^2$, "
          "$m_\\psi^2(z)=m_0^2+g^2f^2x(z)^2$, and the modes of the planar wall in its rest frame obey"))
add(("eq", r"\left[-\frac{d^2}{dz^2}+m_\psi^2(z)\right]u=(\omega^2-k_\parallel^2)\,u,\qquad T=\frac{k_+}{k_-}\frac{1}{|A|^2},\qquad R=\left|\frac{B}{A}\right|^2"))
add(("p", "with a transmitted wave $e^{ik_+z}$ on the false-vacuum side, integrated back across the wall and decomposed into $Ae^{ik_-z}+Be^{-ik_-z}$ on the "
          "true-vacuum side (interface_transport.py). The same equation describes the wave dark matter of Section 3.3, whose mass changes by the factor "
          "$A(x_t)/A(x_h)$ across the wall. The solver is validated on the symmetric kink, for which $m^2=g^2f^2-g^2f^2\\,{\\rm sech}^2(z/2)$ is a "
          "Poeschl--Teller well with the exact transmission $T=\\sinh^2(2\\pi q)/[\\sinh^2(2\\pi q)+\\cos^2(\\frac\\pi2\\sqrt{1+16g^2/\\lambda})]$; the numerical "
          "values agree with it to %s, including the reflectionless coupling $g^2/\\lambda=1/2$, and $T+R=1$ holds to %s (Fig. 4b)."
          % (sci(TJ["max_PT_error"], 0), sci(TJ["max_unitarity_error"], 0))))
add(("p", "Through the wall at the trigger the interior is heavier for every coupling ($m_{\\rm in}/m_{\\rm out}=%.2f$, since $|x_t|>|x_h|$). Modes below the "
          "interior threshold are reflected, and above it $T$ exceeds 0.5 within $0.1\\sqrt\\lambda f$ of the threshold for every coupling tested (Fig. 4a). "
          "The interface is selective, but kinematically: it filters by mass, as any step does, with little additional reflection from its shape. Photons, "
          "neutrinos and baryons do not couple to $\\varphi$ and cross unchanged. A thermal bath of test scalars exerts on a static wall a reflected momentum "
          "flux that is at most %s of the driving pressure at temperatures up to the wall scale, and in the ultrarelativistic limit at most "
          "$\\Delta m^2T^2/24$ per degree of freedom \\cite{BM09}; at the CMB temperature of the trigger epoch this is below %s of the driving pressure."
          % (g1["m_in"] / g1["m_out"], sci(max(r["P_static_over_dV"] for r in TJ["pressure"] if r["T"] <= 0.2), 0), sci(0.31 * TJ["T_cmb_at_trigger_units"] ** 2, 0))))
add(("p", "For the coupled dark matter the mass ratio across the wall is $A(x_t)/A(x_h)=%.2f$, far larger than its kinetic energy, so a slow wall would "
          "reflect all of it, as at the filtering transitions of \\cite{BKL20}. The wall is not slow. Its driving pressure exceeds the largest pressure the "
          "dark matter can exert, $\\rho_*\\Delta A$ (the Boedeker--Moore limit for quanta that are transmitted at high Lorentz factor \\cite{BM09}), by a "
          "factor %.0f, so it runs away; dark matter at rest crosses once the wall's Lorentz factor exceeds $A(x_t)/A(x_h)$, which happens within %.2f mm of the "
          "nucleation point. The fraction of the dark matter swept by the wall before it becomes transparent is $(\\gamma_{\\rm th}R_0/R_{\\rm coll})^3=%s$. "
          "The interface is therefore transparent to the dark matter for all practical purposes, and the transitioned region inherits the dark matter of its "
          "surroundings with the heavier mass of the true vacuum."
          % (DM["A_true"] / DM["A_false"], 1.0 / DM["P_dm_max_over_dV_vac"], DM["R_threshold_mm"], texsci(DM["filtered_fraction"], 0))))
rA = DM["A_true"] / DM["A_false"]
dep_ad, dep_su = DM["A_true"] - DM["A_false"], DM["A_false"] * (1 + rA * rA) / 2 - DM["A_false"]
gam_sud = TW["thickness_mm"] * GT["m_phi_meV"] * 10 / (HBARC_MM := 1.97327e-1)
add(("p", "For wave dark matter the crossing has a second aspect. A coherent field whose mass jumps from $m$ to $rm$, with $r=A(x_t)/A(x_h)=%.2f$, "
          "keeps its number of quanta if the jump is slow compared with its oscillation (energy density $\\times r$) and keeps its amplitude and velocity if "
          "the jump is sudden (energy density $\\times(1+r^2)/2$, and $(1+r^2)/2r=%.2f$ times as many quanta). The wall crosses a dark-matter field of mass "
          "$m_{\\rm dm}$ in a time (thickness)/$\\gamma$, so the crossing becomes sudden once $\\gamma$ exceeds the wall thickness in units of $1/m_{\\rm dm}$, "
          "about %.0f for $m_{\\rm dm}=10m_\\varphi$, which the runaway wall reaches about a centimetre from the nucleation point. Almost all the dark "
          "matter is therefore crossed suddenly, and the energy the wall deposits in it is $%.2f\\rho_*$ instead of $%.2f\\rho_*$; relative to the driving "
          "pressure this is %.2f rather than %.2f, so the wall still runs away. The transitioned region then holds %.0f%% more dark-matter energy than the "
          "adiabatic estimate used in Section 7, an uncertainty of the same kind as the equation of state of the collision products."
          % (rA, (1 + rA * rA) / (2 * rA), gam_sud, dep_su, dep_ad, DM["P_dm_max_over_dV_vac"] * dep_su / dep_ad, DM["P_dm_max_over_dV_vac"],
             100 * ((1 + rA * rA) / (2 * rA) - 1))))
add(("fig", "fig4_selective_transmission.png", "Fig. 4. (a) Transmission of a test scalar through the wall at the trigger against momentum above the "
            "interior threshold, for four couplings. (b) Validation on the symmetric kink: numerical transmission (dots) against the exact Poeschl--Teller "
            "result (lines), including the reflectionless coupling $g^2/\\lambda=1/2$.", 6.3))

# ============================================================================================ 6
add(("h1", "6. Gravitational wall and collision dynamics"))
add(("p", "For a thin wall of tension $\\sigma$ and areal radius $R$ between an interior of density $\\rho_-$ and an exterior of density $\\rho_+$, the Israel "
          "junction condition can be written with the Misner--Sharp mass on each side \\cite{Israel66,BGG}:"))
add(("eq", r"\Gamma_--\Gamma_+=\frac{\sigma R}{2M_{\rm Pl}^2},\qquad \Gamma_\pm=\sqrt{1+\dot R^2-\frac{\rho_\pm R^2}{3M_{\rm Pl}^2}}\ \Rightarrow\ \Gamma_-=\frac{(\rho_+-\rho_-)R}{3\sigma}+\frac{\sigma R}{4M_{\rm Pl}^2}"))
add(("p", "(bubble_gravity.py). From $R_0=3\\sigma/\\Delta V$ (%.2f mm) to the Hubble radius the wall follows the flat-space runaway $\\Gamma=R/R_0$ to "
          "round-off and reaches $\\Gamma=$%s; its self-gravity, $\\sigma R/2\\Gamma$, is %s; and the Misner--Sharp mass it carries equals the vacuum-energy "
          "deficit of the region it encloses (Fig. 5a). The released vacuum energy is therefore carried by the wall, not deposited in the interior, until walls "
          "collide or stall at the boundary of a denser region."
          % (GJ["wall"]["R0_mm"], sci(GJ["wall"]["Gamma_at_RH"], 1), sci(GJ["wall"]["max_grav"], 0))))
add(("p", "The completion of the transition in a connected region is computed from its rate history (bubble_handoff.py): the false-vacuum fraction "
          "$f=e^{-I}$, $I=\\frac{4\\pi}{3}\\int\\Gamma a^3r^3dt$ with $r$ the comoving radius of a bubble nucleated at $t$ \\cite{GW,TWW}, falls to one half at "
          "$a=%.3f$ and to 1%% at $a=%.3f$ (Fig. 6a). The rate grows at $\\beta=d\\ln\\Gamma/dt=%.0f\\,H$, so the transition is completed in %.3f Hubble times and "
          "the bubbles collide at a mean radius $R_{\\rm coll}=(8\\pi)^{1/3}/\\beta=%.4f/H$, comparable to the size of large voids. At collision the walls "
          "have Lorentz factors of order %s."
          % (PC["a_half"], PC["a_one_percent"], PC["beta_over_H"], PC["completion_duration_over_hubble_time"], PC["R_coll_times_H"],
             sci(WB["gamma_at_collision"], 0))))
gmin = min(r["gamma_terminal"] for r in WB["gauge_nlo_terminal_gamma"] if r["T_meV"] < 1.0 and r["g"] <= 1.0)
add(("p", "The junction treatment omits several effects, which are bounded instead (Table 3). The Lorentz-contracted thickness of the wall is %s of its "
          "radius, so the thin-wall description is excellent once the wall has grown. The pressure of the dark matter is at most %.2f of the drive. A dark gauge "
          "boson acquiring mass from $\\varphi$ would add a pressure that grows with the Lorentz factor \\cite{BM17} and stop the runaway; for couplings "
          "$g\\le1$ and any bath colder than 1 meV the terminal Lorentz factor is at least %.0f, still above the threshold %.1f at which the dark matter "
          "crosses, so the transparency of Section 5 is unchanged. Gravitational waves from the collisions carry at most a fraction $\\sim%s$ of the wall "
          "energy, suppressed by $(H/\\beta)^2$, and since $R_{\\rm coll}H\\approx%.3f$ no region of Hubble size is overdense enough to form black holes."
          % (sci(WB["contracted_thickness_over_radius"], 0), WB["dm_pressure_over_vacuum_drive"], gmin, WB["gamma_threshold_dm"],
             texsci(WB["gw_fraction_bound"], 0), WB["pbh_scale_ratio"])))
add(("p", "None of these signals is observable today. The trigger lies in the future, at $a\\approx%.1f$ for the baseline, and a protected past means "
          "precisely that no bubble has nucleated within our past light cone; the imprints of bubble collisions on the microwave background that have been "
          "studied for transitions in the early universe \\cite{Kleban,AJ} therefore do not arise. What is observable now is indirect: the dark energy is "
          "exactly $w=-1$ while the vacuum is held, and the dark matter must be of the kind specified in Section 3."
          % WJ["a_trigger"]))
add(("table", "Table 3. Approximations of the wall dynamics and of the hand-off to the enclosed region, and their bounds (bubble_handoff.py, granularity.py).",
     ["Effect", "Bound or result", "Consequence"],
     [["mean-field source", "wave dark matter, %.0f meV $\\ll m_{\\rm dm}\\ll$ %.0f meV" % tuple(GRJ["wave_window_meV"]), "required (Section 3)"],
      ["completion of the transition", "%.3f Hubble times; $R_{\\rm coll}H=%.4f$" % (PC["completion_duration_over_hubble_time"], PC["R_coll_times_H"]),
       "void by void"],
      ["equation of state of collision products", "$w=1/3$ and $w=0$ bracketed", "turnaround moves by %.0f%%" %
       (100 * abs(W0["a_max_over_a_nuc"][1] / W13["a_max_over_a_nuc"][1] - 1))],
      ["thin wall after contraction", "thickness/radius %s" % texsci(WB["contracted_thickness_over_radius"], 0).join("$$"), "valid"],
      ["self-gravity of the wall", "$\\sigma R/2\\Gamma=%s$" % texsci(GJ["wall"]["max_grav"], 0), "negligible"],
      ["dark-matter pressure", "%.2f of the drive (adiabatic), %.2f (sudden crossing)" % (WB["dm_pressure_over_vacuum_drive"],
                                                                                      WB["dm_pressure_over_vacuum_drive"] * dep_su / dep_ad), "runaway"],
      ["dark-matter energy after crossing", "+%.0f%% if sudden" % (100 * ((1 + rA * rA) / (2 * rA) - 1)), "bracketed"],
      ["gauge-boson (NLO) pressure", "terminal $\\Gamma\\ge%.0f$ ($g\\le1$, $T<1$ meV)" % gmin, "dark matter still crosses"],
      ["gravitational waves", "fraction $\\lesssim%s$" % texsci(WB["gw_fraction_bound"], 0), "negligible"],
      ["black holes from collisions", "$R_{\\rm coll}H=%.4f$" % WB["pbh_scale_ratio"], "none"],
      ["isolated Hubble-scale bubble", "trapped at $t=%.2f$, before the interior bounce at %.2f" % (IS["t_trapped"], IS["t_bounce"]), "black hole"]], None))
add(("fig", "fig5_daughter_evolution.png", "Fig. 5. (a) Lorentz factor of the wall from the Israel junction against the flat-space runaway, and its "
            "self-gravity (right axis). (b) Scale factor of a transitioned region continued formally through an effective bounce (the conditional "
            "continuation of Section 7.3, not supported by the dark matter of Section 3). (c) The same equations for the open region enclosed by an isolated "
            "Hubble-scale bubble, which becomes trapped before its restoration (Section 7.2).", 6.6))

# ============================================================================================ 7
add(("h1", "7. Energy accounting and the fate of the transitioned region"))
add(("h2", "7.1 The region enclosed by the wall"))
nB, tB, rB, bB = ev(B_P, "nucleation"), ev(B_P, "turnaround"), ev(B_P, "restoration"), ev(B_P, "bounce")
tI = ev(B_I, "turnaround")
add(("p", "A wall says nothing by itself about the cosmology inside it; the interior has to be solved. It is evolved as a Friedmann region with radiation, "
          "baryons, dark matter of field-dependent mass $A(x)$ and the scalar in one of its minima: the scalar mass, of order meV, exceeds $H$ by thirty orders "
          "of magnitude, so the field follows its minimum adiabatically, and the stress-energy is $\\rho=\\rho_r+\\rho_b+\\rho_*+U$, "
          "$\\rho+P=\\frac43\\rho_r+\\rho_b+\\rho_*A(x)$, with $U(x,\\theta)$ the vacuum energy including the coupling energy $\\rho_*(A-1)$. The junction fixes "
          "the initial data: the fluid inside keeps the Hubble rate of the parent at nucleation while its density drops by $\\Delta V_{\\rm eff}$, so an isolated "
          "Hubble-scale bubble encloses an open region with $K/a^2=\\Delta V_{\\rm eff}/3$, and a completed transition leaves a flat region whose wall energy is "
          "released at the collisions. Events (nucleation, turnaround, restoration) are found as roots in $\\ln a$ with $H$ from the constraint, so the "
          "constraint holds exactly between events; an independent integration as an initial-value problem reproduces the turnaround and the restoration to "
          "%s in time with a constraint residual of %s."
          % (sci(abs(RP["restoration_time_ivp"] - RP["restoration_time_root"]) / RP["restoration_time_root"], 0), sci(RP["constraint_residual_ivp"], 0))))
add(("p", "In the environmental picture of Section 4.4 the transitioned volume is a network of voids and, later, sheets and filaments, which together fill "
          "most of space, while the haloes keep the held vacuum. Its volume-averaged vacuum energy is negative once it fills most of space, and a region in "
          "which the transition is complete holds the negative true vacuum, turns around at %.2f times the scale factor at nucleation and contracts. When the "
          "dark-matter density in the contraction reaches $\\theta_{\\rm high}$ the true minimum disappears and the field returns to the matter-held minimum, "
          "releasing $%.2f\\,V_f$; the matter term restores the held vacuum, a direct consequence of the density dependence. This restoration needs "
          "dissipation: undamped, the oscillation about the held minimum is an adiabatic invariant that returns to the same amplitude when the density falls "
          "back through $\\theta_{\\rm high}$. The inhomogeneous evolution of a universe with islands of the held vacuum is not computed here; the "
          "volume-averaged statement is that it recollapses."
          % (tB["a_max_over_a_nuc"], rB["dE_released"] / VC)))
add(("h2", "7.2 Isolated Hubble-scale bubbles are cut off"))
add(("p", "A bubble that remains isolated on the Hubble scale, which happens when the transition does not complete, encloses an open region that also "
          "turns around (at %.2f times the nucleation scale factor), but its continuation is not part of the parent universe. Following the wall as a null ray "
          "in the open interior and comparing the Misner--Sharp mass inside it with that of the exterior false-vacuum universe at the same areal radius, the "
          "interior becomes a trapped region ($\\rho_{\\rm in}R^2/3>1$) at $t=%.2f/H_0$ and denser than the exterior at $t=%.2f/H_0$, both before its "
          "restoration ($%.2f$) and before any bounce ($%.2f$) (Fig. 6b). From then on the junction requires $\\Gamma_+$ to change sign: the recollapsing region "
          "is cut off behind a horizon and is seen from outside as a black hole, as for the false-vacuum bubbles of \\cite{BGG}, and whatever happens inside it "
          "is causally disconnected from the parent \\cite{FG,FGG}. Void-sized regions, much smaller than the Hubble radius, are not trapped in this way: their "
          "trapping parameter is smaller by $(R_{\\rm void}H)^2\\sim10^{-4}$, and they follow the averaged evolution."
          % (tI["a_max_over_a_nuc"], IS["t_trapped"], IS["t_interior_denser"], IS["t_restoration"], IS["t_bounce"])))
add(("fig", "fig6_handoff.png", "Fig. 6. (a) Fraction of space in the false vacuum along the dilution history of a connected region, with the completion "
            "criterion $\\Gamma/H^4=9/4\\pi$. (b) Isolated Hubble-scale bubble: density inside the wall and in the exterior, and the trapping parameter "
            "$\\rho_{\\rm in}R^2/3$, with the turnaround, restoration and bounce of the interior.", 6.4))
add(("h2", "7.3 Energy and entropy closure, and why a continuation is not supported"))
add(("p", "Every species crosses the wall, so the comoving entropy of the region is inherited, and it is raised when the wall energy is thermalized at "
          "the collisions (by a factor of about %.0f) and when the vacuum is restored (bubble_entropy.py); it never decreases. The interior bookkeeping "
          "$dE_{\\rm in}=F_E\\,dV+W_{\\rm wall}$, with the swept-in content $F_E=\\rho_+$ and the work $W_{\\rm wall}=-\\Delta V_{\\rm eff}dV$ handed to the "
          "wall, closes to round-off. The decay removes the cosmological horizon of the false vacuum, whose entropy ($%s$) exceeds that of the radiation in a "
          "Hubble volume ($%s$); this is the usual situation for decays out of de Sitter space, accounted for by the action itself."
          % (EJ["entropy"][B_P]["nucleation_jumps"][1], texsci(EJ["horizon"]["S_dS"], 1), texsci(EJ["horizon"]["S_radiation_Hubble_volume"], 1))))
add(("p", "The collision products are not computed; their equation of state is bracketed: if they behave as radiation the region evolves exactly as with "
          "thermalized wall energy, and if they are massive scalar quanta ($w=0$) the turnaround moves from %.2f to %.2f times the nucleation scale factor. "
          "Whether the recollapsing region could be continued through a nonsingular bounce into a new cycle is a separate question. With an effective "
          "loop-quantum-cosmology bounce \\cite{APS06,AS11} the equations do repeat, self-similarly with growth $\\ln G=%.3f$ per cycle ($%.3f$ for $w=0$ "
          "products), but only if the bounce returns the dark matter to the state that holds the vacuum in the next cycle. For a thermal relic a hot bounce "
          "does that through freeze-out; for the wave dark matter that Section 3 requires it does not, because a meV-scale boson thermalized at a hot bounce is "
          "a hot relic, not a cold coherent condensate. The dark matter that allows the trigger is therefore not the dark matter that would allow a cyclic "
          "continuation, and the latter is not supported in this branch."
          % (W13["a_max_over_a_nuc"][1], W0["a_max_over_a_nuc"][1], LNG["B"], float(np.mean(W0["lnG"])))))

# ============================================================================================ 8
add(("h1", "8. Parameter-space analysis and numerical controls"))
add(("p", "Controls separate the contribution of each ingredient (bubble_to_cycle.py, granularity.py; Table 4). Without the coupling (D) there is no "
          "density dependence, and a single transition that has neither already happened in our past light cone nor fails to complete requires $B_0$ in a "
          "window of width %.0f around %.0f. With the coupling but a granular source, particle dark matter or wave dark matter lighter than the field, the "
          "rate is the zero-density rate, and the vacuum decays in the past. With the computed transmission (C) the region is identical to that with every "
          "species transmitted (B) to %s: the selectivity of the wall changes nothing, and a hypothetical wall that keeps all dark matter out (C*) differs "
          "only in the conditional continuation."
          % (DW[1] - DW[0], 0.5 * (DW[0] + DW[1]), sci(DM["filtered_fraction"], 0))))
rows4 = [["D", "no matter coupling", "constant rate", "window $%.0f\\le B_0\\le%.0f$" % tuple(DW), "fine tuning"],
         ["G1", "coupling, particle dark matter", "zero-density rate between particles", "$e^{%.0f}$ bubbles in the past" %
          GRJ["particle_dm"]["ln_expected_bubbles_past_light_cone"], "excluded"],
         ["G2", "coupling, wave dark matter $m_{\\rm dm}\\ll m_\\varphi$", "zero-density rate at the nodes", "as G1", "excluded"],
         ["B", "coupling, wave dark matter in the window, all species transmitted", "local mean field", "past clear in voids for $\\bar\\theta_0\\gtrsim%.1f\\,"
          "\\theta_{\\rm deg}$" % VOID[0.1]["min_f0_mean"], "viable"],
         ["C", "as B with the computed transmission", "local mean field", "identical to B", "viable"],
         ["C*", "fully selective wall (hypothetical)", "local mean field", "differs only in a continuation", "--"]]
add(("table", "Table 4. Controls ($\\eta=0.30$, $q=1$, $B_0=%.0f$; voids at %.0f%% of the mean density)." % (WJ["base"]["B0"], 100 * DV_),
     ["Model", "Description", "Release", "Past light cone", "Outcome"], rows4, None))
add(("p", "In the plane of $B_0$ and the local coupling today in the emptiest regions, the region in which the past light cone is clear and the transition "
          "completes extends from $B_0=%.0f$ to $%.0f$ for suitable couplings and covers %.0f%% of the plane scanned (Fig. 7); above the upper end the vacuum "
          "rate itself is too small for completion. A survey over stated ranges (viability_mc.py; %d draws, log-uniform $B_0\\in[50,1000]$, local "
          "$\\theta_0/\\theta_{\\rm deg}\\in[10^{-3},10]$, uniform $q\\in[-1,3]$, $\\eta\\in\\{0.20,0.30,0.36\\}$) finds a clear past light cone in %.0f%% of draws, "
          "completion in %.0f%% and a negative true vacuum, needed for recollapse, in %.0f%%; the three hold together in %.0f%%. These fractions describe the "
          "sample space of the model, not a probability for our universe, and they presuppose the wave dark matter of Section 3; for particle dark matter "
          "the fraction is zero."
          % (CJ["map"]["B0_range_viable"][0], CJ["map"]["B0_range_viable"][1], 100 * CJ["map"]["fraction_viable"], TM["n"],
             100 * TM["single"]["T1 past light cone clear"], 100 * TM["single"]["T2 completes"], 100 * TM["single"]["T3 recollapse"],
             100 * TM["joint_background"])))
add(("p", "All %d checks of interface_validation.py pass (Table S1): the analytic degeneracy and window edges, the kink tension and thickness, the bounce "
          "actions against an independent table, the thin-wall limit, shooting against collocation, the transmission against Poeschl--Teller and the sharp "
          "step, unitarity and tolerance convergence, the flat-space limit of the junction and the energy carried by the wall, the constraint and continuity "
          "of the region evolution, entropy monotonicity, the controls, the completion time and the independent integration. The granularity estimates "
          "(granularity.py) are closed-form and are tabulated in Tables S4 and S5." % VJ["n"]))
add(("fig", "fig7_viability_map.png", "Fig. 7. Outcome over the vacuum action $B_0$ and the local coupling today in the emptiest regions ($\\eta=0.30$, "
            "$q=1$): already decayed in the past light cone, no completed transition, or a completed transition; the red band is the window without "
            "coupling (D).", 4.6))

# ============================================================================================ 9
add(("h1", "9. Implications, limitations and conclusions"))
add(("p", "A coupling of a scalar to dark matter centred on a false vacuum makes the stability of that vacuum a function of the dark-matter density, and "
          "dilution turns the held vacuum from the true into a metastable one at $\\theta_{\\rm deg}=\\eta/2$. For a vacuum at the dark-energy scale this "
          "mean-field statement holds only if the dark matter looks smooth to the field on a scale below a tenth of a millimetre and a time below a picosecond. "
          "Particle dark matter never does, and with it the vacuum would have decayed long ago; the coincidence of the dark-matter and dark-energy densities, "
          "which motivates a dilution trigger, is what forbids it. Wave dark matter heavier than the field and still classical, roughly 3 meV to 0.25 eV, "
          "does, and then the decay is environmental. It begins in the emptiest voids, the past is protected only if the mean coupling today is close to the "
          "degeneracy density, the transition advances from voids into sheets and filaments as the universe expands, and galactic haloes keep the held "
          "vacuum. Within the voids the wall is a fraction of a millimetre thick, runs away, becomes transparent to the dark matter and carries the released "
          "energy; transitioned regions recollapse in their negative vacuum, and isolated Hubble-scale bubbles end as black holes."))
add(("p", "The observational implications follow from the environment. The dark energy measured inside haloes is unchanged until the global recollapse, "
          "so local tests see none of the transition, while the volume-averaged expansion responds to the voids. The dark matter must be a nonthermal boson in "
          "the meV to sub-eV range with a coupling to the scalar of strength $M\\approx f$, which changes its mass by a factor %.1f between the vacua; such a "
          "coupling acts only over the Compton length of the scalar and is invisible to astrophysical fifth-force tests. A detection of thermal-relic dark "
          "matter, or of dark matter heavier than about an eV, would exclude the density-triggered decay of a dark-energy-scale vacuum in this class of models."
          % (DM["A_true"] / DM["A_false"])))
add(("p", "The limitations are specific. The local density is treated as constant along the past light cone of a region, which is conservative for voids; "
          "the void statistics enter only through the depth of the emptiest regions; the advance of the walls into denser regions and the islands of held "
          "vacuum are described through the degeneracy contour, not simulated; and the backreaction of the scalar on the wave dark matter (a parametric "
          "coupling between two oscillating fields) is not computed. The one-loop prefactor is estimated, not computed, and enters the completion time only "
          "logarithmically; the collision products are bracketed by their equation of state; and the next-to-leading friction is estimated for a hypothetical "
          "dark gauge boson only."))
add(("p", "The results do not require a cosmological bounce, and they do not support one in this branch. They establish when a density-dependent vacuum "
          "can decay at all, where it decays first, what its walls transmit and how they move, and what the regions they enclose become. The homogeneous, "
          "classical version of the trigger, which does not face the granularity bound because its field has a cosmological Compton length, is examined in the "
          "companion paper \\cite{Companion}."))

# ============================================================================================ declarations
add(("h1", "Declarations"))
add(("pb", "Funding.", "This research received no specific grant from funding agencies in the public, commercial, or not-for-profit sectors."))
add(("pb", "Declaration of competing interest.", "The author declares no known competing financial interests or personal relationships that could have "
        "influenced the work reported in this paper."))
add(("pb", "CRediT authorship contribution statement.", "Leon Sandler: Conceptualization, Methodology, Software, Formal analysis, Investigation, Writing -- "
        "original draft, Writing -- review and editing."))
add(("pb", "Data and code availability.", "All code (Python, NumPy, SciPy, Matplotlib \\cite{NumPy,SciPy}) and the data behind every figure and table are at %s "
        "(MIT license) and archived at https://doi.org/%s; this manuscript is archived at https://doi.org/%s. Input tables shared with the companion paper "
        "are documented as such in the README." % (GITHUB, DOIS.get("code", "(to be assigned)"), DOIS.get("paper", "(to be assigned)"))))
add(("pb", "Related manuscript.", "A companion manuscript on the homogeneous version of the trigger and its continuation into cycles is submitted to The "
        "European Physical Journal C%s; it shares the class of models but none of the results reported here."
        % ((", archived at https://doi.org/%s" % COMP["paper"]) if COMP.get("paper") else "")))
add(("pb", "Declaration of generative AI and AI-assisted technologies in the writing process.", "During the preparation of this work the author used Claude "
        "(Anthropic) to assist with the numerical analysis, code, literature checking and drafting. After using this tool, the author reviewed and edited the "
        "content as needed and takes full responsibility for the content of the publication."))
add(("h1", "References"))

order = []
for blk in B:
    chunks = [blk[-1]] if blk[0] in ("p", "abstract", "pb") else (list(blk[1]) if blk[0] == "bullets" else
                                                                  ([str(c) for row in blk[3] for c in row] if blk[0] == "table" else []))
    for ch in chunks:
        for m in re.finditer(r"\\cite\{([^}]+)\}", ch):
            for key in m.group(1).split(","):
                if key not in order:
                    order.append(key)
undefined = [k for k in order if k not in REFS]
if undefined:
    sys.exit("undefined citation keys %s" % undefined)
numbering = {k: i + 1 for i, k in enumerate(order)}
rep = lambda s: re.sub(r"\\cite\{([^}]+)\}", lambda m: "[" + ", ".join(str(numbering[x]) for x in m.group(1).split(",")) + "]", s)
for k in order:
    add(("ref", "[%d] %s" % (numbering[k], REFS[k])))
B2 = []
for b in B:
    if b[0] in ("p", "abstract", "note"):
        B2.append((b[0], rep(b[1])))
    elif b[0] == "pb":
        B2.append(("pb", b[1], rep(b[2])))
    elif b[0] == "table":
        B2.append(("table", b[1], b[2], [[rep(str(c)) for c in row] for row in b[3]], b[4]))
    else:
        B2.append(b)
B = [CE.block(b) for b in B2]
BASENAME = "PDU_Vacuum_Bubbles"


def supplement():
    S = []
    a = S.append
    a(("title", "Supplementary material: " + TITLE))
    a(("author", "Leon Sandler\nIndependent Researcher, Northbrook, Illinois, USA"))
    a(("h1", "S1. Methods"))
    a(("p", "granularity.py evaluates the scales of Section 3 (the field mass at the held vacuum, the Compton length, the dark-matter density, the particle, "
            "classicality and coherence bounds), the zero-density rate and past-light-cone integral seen between particles, the structural ratio over the "
            "bounce tables and toward $\\eta_c$ (with new bounces computed by vacuum_wall.bounce_profile), the protection of the past in regions of reduced "
            "density, and the void-by-void completion delay. vacuum_phase_map.py maps the stationary points over density and tilt and evaluates the coherence "
            "of the symmetron release. vacuum_wall.py solves the O(4) bounce with its profile, the thin-wall tension with the energy difference removed, the "
            "static wall at the degeneracy density, and the nucleation history with the actual Hubble rate. interface_transport.py integrates the mode equation "
            "across the planar wall with complex amplitudes (DOP853, rtol $10^{-10}$) and computes the reflected momentum flux of a Bose bath with the "
            "transverse integral in closed form. bubble_gravity.py evaluates the Israel junction and evolves the enclosed regions by root finding in $\\ln a$; "
            "bubble_to_cycle.py runs the controls and the viability map; bubble_entropy.py the energy and entropy closure; bubble_handoff.py the percolation "
            "integral, the bounds of Table 3, the isolated-bubble test, the collision-product variants and the independent integration; viability_mc.py the "
            "survey; interface_validation.py the checks of Table S1."))
    a(("table", "Table S1. Validation (interface_validation.py).", ["Check", "Result"], IB.validation_rows(), None))
    rows = [["%.2f" % r["gh"], "%.2f" % r["T"], texsci(max(r["P_static_over_dV"], 1e-300), 1).join("$$"), texsci(max(r["P_BM_over_dV"], 1e-300), 1).join("$$")]
            for r in TJ["pressure"]]
    a(("table", "Table S2. Pressure of a thermal test-scalar bath on the wall at the trigger relative to the driving pressure: static wall and "
                "ultrarelativistic limit $\\Delta m^2T^2/24$. Temperatures in units of $\\sqrt\\lambda f$; the CMB at the trigger is %.3f."
                % TJ["T_cmb_at_trigger_units"], ["$g/\\sqrt\\lambda$", "T", "static", "relativistic limit"], rows, None))
    rows = [[v["variant"], ", ".join("%.4f" % x for x in v["lnG"]), "%.2f" % v["T_b_MeV"], ", ".join("%.3f" % x for x in v["a_max_over_a_nuc"])]
            for v in HJ["variants"]]
    a(("table", "Table S3. Formal continuation variants (not supported by the dark matter of Section 3; Section 7.3): growth per cycle, bounce temperature "
                "($\\rho_r=T^4$) and turnaround relative to nucleation.", ["Variant", "$\\ln G$", "$T_b$ (MeV)", "$a_{\\max}/a_{\\rm nuc}$"], rows, None))
    rows = [["%.3f" % r["eta"], "%.1f" % r["theta_over_deg"], "%.1f" % r["q"], "%.1f" % r["b"], texsci(r["barrier"], 1).join("$$"), "%.3f" % r["lam"],
             "%.2f" % r["ratio"]] for r in GRJ["near_eta_c"]]
    a(("table", "Table S4. Approach to $\\eta_c$: bounce action, barrier, $\\lambda=b/B$ and $\\rho_{\\rm dm}/m_\\varphi^4$ at the trigger epoch.",
       ["$\\eta$", "$\\theta/\\theta_{\\rm deg}$", "$q$", "$b$", "barrier", "$\\lambda$", "$\\rho_{\\rm dm}/m_\\varphi^4$"], rows, None))
    rows = [["%.1f" % v["one_plus_delta"], "%.1f" % v["baseline_lnI"], "yes" if v["baseline_past_clear"] else "no", "%.3f" % v["min_f0_mean"]]
            for v in GRJ["voids"]]
    a(("table", "Table S5. Regions of reduced density: ln of the expected number of bubbles in the past light cone for a mean coupling "
                "$\\theta_0=%.1f\\,\\theta_{\\rm deg}$, and the smallest mean coupling for which the past is clear." % F0_LOC,
       ["$1+\\delta$", "ln expected bubbles", "past clear", "minimum mean $\\theta_0/\\theta_{\\rm deg}$"], rows, None))
    a(("h1", "S2. Reproduction"))
    a(("p", "Environment of the published run: Python %s, NumPy %s, SciPy %s (Windows). From the code directory: vacuum_phase_map.py, vacuum_wall.py, "
            "interface_transport.py (about 25 minutes), bubble_gravity.py, bubble_to_cycle.py, bubble_handoff.py, viability_mc.py (about 15 minutes), "
            "bubble_entropy.py, interface_validation.py, granularity.py, make_int_figures.py, make_gran_figure.py, paper2_build.py. Input tables shared with "
            "the companion paper (bounce_tables.json, thermal_tables.json, survival.json, density_bounce.json, coupling.json, continuous_cycles.json, "
            "cycle_map.json) are included as data and documented in the README. Code: %s." % (platform.python_version(), np.__version__, scipy.__version__,
                                                                                           GITHUB)))
    S = [CE.block(b) for b in S]
    RN.render_docx(S, os.path.join(OUT, BASENAME + "_Supplement.docx"), FIGDIR)


def letter():
    d = Document()
    st = d.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(11); st.paragraph_format.space_after = Pt(8)
    for s in d.sections:
        s.left_margin = s.right_margin = Inches(1.1); s.top_margin = s.bottom_margin = Inches(1.0)

    def para(text, bold=False, after=None):
        p = d.add_paragraph(); r = p.add_run(CE.ce(text)); r.bold = bold
        if after is not None:
            p.paragraph_format.space_after = Pt(after)
    para("Leon Sandler", True, 0)
    for t in ("Independent Researcher", "Northbrook, Illinois, USA", "sandler.leon@gmail.com", "ORCID 0009-0007-4584-808X"):
        para(t, after=0)
    para("")
    para(datetime.date.today().strftime("%d %B %Y").lstrip("0"))
    para("The Editors", after=0); para("Physics of the Dark Universe", after=0); para("Elsevier")
    para("Dear Editors,")
    para("I submit the manuscript “%s” as a research article for consideration in Physics of the Dark Universe." % TITLE)
    para("The paper asks whether a declining dark-matter density can switch on the decay of a vacuum, and what the resulting vacuum walls do. For a scalar "
         "coupled to dark matter it maps where the held vacuum is stable and where it becomes metastable, and then asks a question that the literature on "
         "environment-dependent vacua has not had to face: when does the field feel the smooth mean density that such an argument assumes? For a vacuum at "
         "the dark-energy scale the field responds within about 0.07 mm, inside which the dark-matter energy is only a few percent of the field mass to the "
         "fourth power. Particle dark matter of any mass is then a set of point sources between which the vacuum decays at its zero-density rate, so a "
         "density-triggered decay is impossible with it; the coincidence of the dark-matter and dark-energy densities, which motivates such a trigger, is "
         "what forbids it. The mean field is restored only for wave dark matter between about 3 meV and 0.25 eV.")
    para("For that dark matter the transition is environmental: it begins in the emptiest voids, the past is protected only if the mean coupling today is "
         "close to the degeneracy density, and galactic haloes keep the old vacuum. The paper then computes the Coleman-De Luccia bounce and the wall, the "
         "transmission of fields and of the dark matter through the wall (which runs away within a millimetre and filters a fraction of order 1e-88), the "
         "Israel junction, the completion of the transition and the fate of the enclosed regions; an isolated Hubble-scale bubble ends as a black hole. A "
         "continuation through a bounce is shown not to be supported, because it would need a thermal relic.")
    para("Related submission. The homogeneous, classical version of the trigger and its continuation into cycles are the subject of a separate manuscript "
         "submitted to The European Physical Journal C%s. Both grew out of one earlier, broader draft that a review recommended dividing; the granularity "
         "bound reported here does not affect that paper, whose field has a cosmological Compton length. The two share the class of models but not their "
         "results, and I will provide the companion manuscript on request."
         % ((" (preprint https://doi.org/%s)" % COMP["paper"]) if COMP.get("paper") else ""))
    para("The analysis is fully reproducible: the code (Python, MIT license) and the data behind every figure and table are at %s and archived at "
         "https://doi.org/%s; the manuscript is archived at https://doi.org/%s, and it is built only if the %d validation checks pass."
         % (GITHUB, DOIS.get("code", "(to be assigned)"), DOIS.get("paper", "(to be assigned)"), VJ["n"]))
    para("A generative AI tool (Claude, Anthropic) assisted with the numerical analysis, code, literature checking and drafting; I reviewed all of it and take "
         "full responsibility for the content, as declared in the manuscript. The manuscript is original, has not been published elsewhere and is not under "
         "consideration at any other journal. I am the sole author, I have no competing interests, and the work received no external funding.")
    para("Thank you for considering this submission.")
    para("Sincerely,")
    para("Leon Sandler", True)
    d.save(os.path.join(OUT, "PDU_Cover_Letter.docx"))


if __name__ == "__main__":
    for h in HIGHLIGHTS:
        assert len(h) <= 85, (len(h), h)
    os.makedirs(OUT, exist_ok=True)
    Bd = B[:3] + [("h1", "Highlights"), ("bullets", HIGHLIGHTS)] + B[3:]
    RN.render_docx(Bd, os.path.join(OUT, BASENAME + ".docx"), FIGDIR)
    RN.render_tex_els(B, os.path.join(OUT, BASENAME + ".tex"), "Physics of the Dark Universe", HIGHLIGHTS, "sandler.leon@gmail.com",
                      "organization={Independent Researcher},city={Northbrook},state={Illinois},country={USA}")
    supplement()
    letter()
    words = sum(len(b[1].split()) for b in B if b[0] in ("p", "abstract")) + sum(len(x.split()) for b in B if b[0] == "bullets" for x in b[1])
    print("Paper 2: %d words, %d figures, %d tables, %d equations, %d references" % (words, sum(1 for b in B if b[0] == "fig"),
                                                                                    sum(1 for b in B if b[0] == "table"), sum(1 for b in B if b[0] == "eq"), len(order)))
