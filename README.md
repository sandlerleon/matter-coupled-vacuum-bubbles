# Density-Dependent Vacuum Bubble Nucleation in Matter-Coupled Dark-Sector Cosmology

**Dark-Matter Granularity, Wall Dynamics and Transmission**

Leon Sandler, Independent Researcher — sandler.leon@gmail.com — ORCID 0009-0007-4584-808X

Code, data and manuscript for a paper prepared for **Physics of the Dark Universe** (Elsevier). Manuscript:
`manuscript/PDU_Vacuum_Bubbles.{docx,tex,pdf}` (elsarticle), supplement `manuscript/PDU_Vacuum_Bubbles_Supplement.docx`, cover letter
`manuscript/PDU_Cover_Letter.docx`. Archived on Zenodo: paper 10.5281/zenodo.23850096, code 10.5281/zenodo.23850099 (version 1.1.1; all versions:
10.5281/zenodo.23841873 / 10.5281/zenodo.23841877).

## What the paper shows

A scalar coupled to dark matter through A(φ) = 1 + (φ − f)²/2M², centred on a false vacuum, is held there by the matter term while the dark-matter
density is high. In the mean-field description, dilution makes the held vacuum metastable below θ_deg = η/2 and Coleman–De Luccia tunneling switches on.

* **Granularity.** For a vacuum at the dark-energy scale the field's Compton length is ≈ 0.07 mm and the dark-matter energy in a Compton volume is
  ≈ 0.05 m_φ⁴. Spatially resolved tests — the static field in a cell around one source, and the decay rate with Poisson-distributed sources at first
  order in each source's coupling — show that with particle dark matter the vacuum decays at its zero-density rate (for 1 MeV particles ~10⁻⁸ of the
  mean-field suppression survives). The mean field holds for **wave dark matter heavier than the field and still classical, ≈ 3 meV – 0.25 eV**.
* **Wave dark matter and the wall (1+1 D, backreaction).** The wall's Lorentz factor is 7.7% below the mean-field run; the dark matter receives the
  sudden-crossing energy (4.07 vs (1+r²)/2 = 4.19 per swept density) and trails the wall at 0.75c (predicted (r²−1)/(r²+1) = 0.76).
* **Environment.** Decay begins in the emptiest voids; a region completes its transition at a ∝ (1+δ)^{1/3}; galactic halos keep the held vacuum.
* **Fate.** Transitioned regions recollapse; an isolated Hubble-scale bubble develops a trapped surface (apparent horizon) before its interior bounces,
  insensitive to the exterior cosmology. A bounce continuation would need a thermal relic and is not supported.
* **Optional neutrino coupling.** Invisible to DUNE/Hyper-K (constant shift, ΔP = 0); environment- and redshift-dependent neutrino mass of order
  0.3–0.6 y meV; relic neutrinos do not disturb the trigger for y ≲ 30.
* **Surveys.** Background conditions 20% (broad priors) / 41% (focused); with the mean-field and runaway conditions 0.5% / 12%.

## Scripts and results

| Result | Paper | Scripts |
|---|---|---|
| Phase structure, degeneracy θ = η/2, symmetron contrast, coupling check | 2 | `vacuum_phase_map.py` |
| Granularity criterion, wave-dark-matter window, structural ratio and approach to η_c, voids, void-by-void completion | 3, 4.4 | `granularity.py` |
| Static field around discrete sources; granular decay rate (Poisson, first order) | 3.5 | `granular_field.py` |
| Nonperturbative bounce with a source on its axis; self-screening of a source (thin shell) | 3.5 | `nonpert_bounce.py` |
| Granular rate with self-screened particles | 3.5 | `screened_rate.py` |
| Three-dimensional bounce with a row of discrete sources against the same coupling smeared | 3.5 | `multisource_bounce.py` |
| Transition order with evolving void contrasts | 4.4 | `env_evolution.py` |
| O(4) bounce, wall tension and thickness, thin-wall comparison, nucleation history | 4 | `vacuum_wall.py`, `cdl_lib.py`, `thermal.py` |
| Transmission through the wall (Pöschl–Teller validated), thermal pressure, dark-matter runaway | 5 | `interface_transport.py` |
| Wall through wave dark matter in 1+1 D with backreaction | 5 | `wave_wall_1d.py` |
| Israel junction, region evolution, controls, viability map | 6, 7, 8 | `bubble_gravity.py`, `bubble_to_cycle.py` |
| Percolation, bounds, isolated bubble, variants, independent integration | 6, 7 | `bubble_handoff.py` |
| Sudden-crossing energy in the region; trapping sensitivity | 7 | `region_sensitivity.py` |
| Energy and entropy closure | 7.3 | `bubble_entropy.py` |
| Parameter surveys; transition order a ∝ (1+δ)^{1/3} | 8, 4.4 | `viability_mc.py`, `viability_extended.py` |
| Optional neutrino coupling; early-universe scalar mass and CMB free-streaming bound | 9 | `neutrino_coupling.py` |
| 40 validation checks (the manuscript is built only if they pass) | 8, S1 | `interface_validation.py` |

## Reproduction

Python 3.11 with the packages in `requirements.txt`; pdflatex with `elsarticle.cls` for the PDF. From `code/`:

```
python vacuum_phase_map.py && python vacuum_wall.py && python interface_transport.py      # ~25 min for the transport
python bubble_gravity.py && python bubble_to_cycle.py && python bubble_handoff.py && python bubble_entropy.py
python viability_mc.py && python interface_validation.py && python granularity.py
python granular_field.py && python nonpert_bounce.py && python nonpert_bounce.py convergence && python nonpert_bounce.py saturation && python multisource_bounce.py && python screened_rate.py && python wave_wall_1d.py && python region_sensitivity.py && python viability_extended.py && python env_evolution.py && python neutrino_coupling.py
python make_int_figures.py && python make_gran_figure.py && python make_resolved_figures.py && python paper2_build.py
cd ../manuscript && pdflatex PDU_Vacuum_Bubbles.tex && pdflatex PDU_Vacuum_Bubbles.tex
```

**Clean reproducibility check (2026-10-31).** A fresh clone in a new virtual environment (Python 3.11, NumPy 2.4, SciPy 1.17, packages from `requirements.txt`) reran `granular_field.py`, `wave_wall_1d.py`, `region_sensitivity.py`, `viability_extended.py`, `neutrino_coupling.py`, `interface_validation.py`, the figure scripts and `paper2_build.py`: all 10,382 numbers in the regenerated JSON files are identical to the committed ones and the regenerated LaTeX manuscript is identical.

Input tables shared with the companion paper (computed by scripts in `sandlerleon/late-time-vacuum-decay`, tag v1.8.2) are included as data:
`bounce_tables.json`, `thermal_tables.json`, `survival.json`, `density_bounce.json`, `coupling.json`, `continuous_cycles.json`, `cycle_map.json`.

## Companion paper

The homogeneous (symmetron) version of the trigger, whose field has a cosmological Compton length and is therefore not affected by the granularity
criterion, and its continuation into cycles are treated in "Matter-Dilution-Triggered Recollapse and Conditional Nonsingular Cycles in a Metastable
Dark-Energy Model" (General Relativity and Gravitation; code `sandlerleon/late-time-vacuum-decay`). The two papers share the class of models but not
their results.

## License

MIT (code and data); the manuscript text is © the author.
