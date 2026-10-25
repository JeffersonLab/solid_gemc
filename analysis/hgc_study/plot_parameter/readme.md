# plot_parameter

Plots the optical parameters of the Heavy Gas Cherenkov (HGC) and Light Gas
Cherenkov (LGC): PMT quantum efficiency (QE), mirror reflectivity, gas
refractive index, PMT window refractive index and gas absorption length.

## how to run

    root -b -q plot_parameter.C

It writes into the current directory:

    plot_parameter_QE.pdf
    plot_parameter_reflectivity.pdf
    plot_parameter_RefractiveIndex.pdf
    plot_parameter_RefractiveIndex_LGC_UVwin.pdf
    plot_parameter_abs_length.pdf

(pdf instead of png because ROOT in the jlab_ce 2.5 container cannot write png.)

## which curves the simulation uses

Checked against the generated text files gemc reads:
`../../../geometry/hgc_moved/solid_SIDIS_hgc__*_Original.txt` and
`../../../geometry/lgc_moved/lg_cherenkov__*_{SIDIS,PVDIS}.txt`.
The values in the plot match those files exactly.

| plot | curve | HGC | LGC |
|---|---|---|---|
| QE | H12700_03_WLS_meas | used, in analysis only (`../analysis_tree_solid_hgc.C`, `eff_PMT`) | - |
| QE | LGC_H12700 (= WLS_meas rounded, plus 2 UV points) | - | used in gemc: `solid_lgc` hit processor reads `EFFICIENCY` of photocathode `SL_H12700` |
| QE | H12700_03, H12700_03_Hamamatsu | not used | not used |
| QE | QE_H8500_03 | on surface `SL_HGC_mirror_pmt`, which no volume uses | `@efficiency_pmt` in `solid_lgc_mirror.pl`, never assigned |
| reflectivity | reflectivity / LGC mirror, Winston cone (identical) | used: `SL_HGC_mirror` (60 mirrors) | used: `LGC_Mirror`, `LGC_WinstonCone` |
| refractive index | n C4F8O 1.5atm | used: gas `SL_HGC_C4F8O` | - |
| refractive index | n C4F8O 1.0atm, n C4F10 1.5atm, n C4F10 1.0atm | not used | not used |
| refractive index | n LGC CO2 1.0atm | - | used: gas `SL_LGCCgas_SIDIS` (both SIDIS and PVDIS) |
| UV window index | LGC PMT UV window | - | used: `SL_H12700UVwin` |
| absorption length | abs_length / LGC CO2 1.0atm (identical, CO2 table) | none: gas absorption length is commented out in `solid_SIDIS_hgc_materials.pl` | used: `CO2_1atm_AbsLen_alt` |

## what it means for analysis

- HGC: the `solid_hgc` hit processor records raw optical-photon hits
  (pid -22) on the PMT surface, which is made of the gas `SL_HGC_C4F8O`
  (`$material_pmt_surface` in `solid_SIDIS_hgc_geometry.pl`). No QE is applied
  in gemc; photoelectrons come from applying H12700_03_WLS_meas in
  `../analysis.C`.
- LGC: the `nphe` branch of `solid_lgc` already has QE applied in gemc. Do not
  apply QE again in analysis.

## LGC PMT UV window

In `../../../geometry/lgc_moved/geometries/build_pmts.pl` every photon reaching
a pixel crosses gas -> window -> gas: a 1.5 mm `SL_H12700UVwin` window
(n = 1.47-1.52), a 100 um CO2 gap, then the 0.2 um pixel `SL_H12700`
(sensitive, n = CO2 index). The window has no absorption length, so it
absorbs nothing, but with no optical surface defined Geant4 applies Fresnel
reflection at both faces: about 4% per face at normal incidence, more at
oblique angles.

A/B test (2026-10-25): SIDIS moved solenoid + LGC only, hall vacuum,
e- 3 GeV, theta 8-15 deg, phi 0-360 deg, vertex z = -350 cm, 500 events per
variant. B is a copy of `lg_cherenkov__*_SIDIS.txt` with the refractive index
of `SL_H12700UVwin` set to the CO2 index of `SL_LGCCgas_SIDIS` (no reflection).

| variant | mean npe/event |
|---|---|
| A, nominal window | 48.22 +- 0.69 |
| B, window n = CO2 | 53.47 +- 0.74 |

A/B = 0.902 +- 0.018: the window removes 9.8% +- 1.8% of LGC photoelectrons.
The H12700 QE from Hamamatsu is quoted per photon incident on the window, so
it already includes the window's reflection; the nominal LGC likely
double-counts this loss.

## other optical differences between HGC and LGC

Same in both: all mirror and cone surfaces are `dielectric_metal`, `polished`,
`unified`, `SkinSurface` with no roughness (sigma alpha) set, and share one
reflectivity table; the QE values are the same (H12700_03_WLS_meas); neither
gas has Rayleigh scattering.

| | HGC | LGC |
|---|---|---|
| QE applied | in analysis: photon kept with QE x 0.8 packing (`../analysis_tree_solid_hgc.C`) | in gemc by `solid_lgc` (`nphe`), then x 0.8 packing in `../analysis_tree_solid_lgc.C` |
| photodetector | one 21.3 x 21.3 cm sensitive gas plane per sector, 20 um, no window; pixels made in analysis | 3 x 3 PMTs x 64 pixels (6.48 mm, 4 um gaps) per sector, about 15.6 cm square, behind window + 100 um gap |
| behind sensor | `Kryptonite` (gemc kills every track in it) | `LGC_PMT_LightStop` surface, reflectivity 0, 2 cm behind pixels |
| cone | pyramidal `Trd`, 45.8 -> 23.3 cm square, 16 cm long, G4_Fe, mirror surface | straight round `Cons`, 30 cm long, glass, mu-metal shield |
| mirror | 1 sphere per sector, R = 210 cm, CFRP | 2 sphere segments per sector, R = 362.5 and 156 cm, RCFP |
| spectral range | 2.04-6.20 eV (200-607 nm) | 2.04-6.89 eV (180-607 nm) for materials; mirror tables stop at 6.2 eV |

So neither raw tree is a final photoelectron count: `solid_hgc` has no QE and
no packing factor, `solid_lgc` `nphe` has QE but no packing factor.

## LGC 180-200 nm photons

Cherenkov photons are only generated inside the gas refractive-index range, so
LGC makes photons down to 180 nm while HGC stops at 200 nm. Above 6.2 eV the LGC
mirror and cone reflectivity tables have ended and Geant4 uses the last value
(0.76); QE there is 0.331; CO2 absorbs strongly (absorption length 0.75 m at
6.53 eV, 0 at 6.89 eV).

A/B test (2026-10-25), same setup as the UV window test, 3000 events per
variant. B removes the 6.53 and 6.89 eV points from every LGC material table
(`SL_LGCCgas_SIDIS`, `SL_N2gas_SIDIS`, `SL_H12700UVwin`, `SL_H12700`), so all
tables end at 6.2 eV = 200 nm.

| variant | mean npe/event |
|---|---|
| A, tables to 180 nm (nominal) | 48.16 +- 0.27 |
| B, tables cut at 200 nm | 46.56 +- 0.27 |

B/A = 0.967 +- 0.008: photons between 180 and 200 nm give 3.3% +- 0.8% of LGC
photoelectrons. The tables hold no measured UV-glass transmission or mirror
reflectivity below 200 nm, so this 3.3% is not physically grounded.

## open questions

- HGC gas has no absorption length, while LGC gas (CO2) has one.
- HGC PMT surface is plain gas, so the H8500 QE on `SL_HGC_mirror_pmt` has no effect.
- LGC UV window Fresnel loss (about 10%) on top of a QE that already includes
  the window; consider an index-matched window or an optical surface.
- LGC 180-200 nm photons (3.3%): stop LGC tables at 200 nm like HGC, or add
  measured window transmission and mirror reflectivity down to 180 nm.
