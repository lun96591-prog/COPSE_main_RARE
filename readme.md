# COPSE

Carbon-Oxygen-Phosphorus-Sulfur-Evolution

COPSE was developed by Noam Bergman in 2003. Several versions now exist. This version of the COPSE code is lighlty updated from Tostevin and Mills (2020) Interface Focus and is the most recent version of the code available (as far as I know).

COPSE has largely been superceded by the SCION model, which can compute spatial surface processes. See https://github.com/bjwmills/SCION for more details.


## Contents:

COPSE_frontend.m // Model frontend file // run "COPSE_frontend(0)" to solve model and plot results // sets parameters and starting values, runs solver

COPSE_equations.m // Model equations file // do not run this code directly // contains flux and reservoir calculations

COPSE_sens.m // sensitivity analysis // individual sensitivity variables are contained in frontend and equations files // call this script to run sensitivity and plot

COPSE_plot.m // called by frontend to plot results // do not run this code directly

COPSE_plot_sens.m // called by sensitivity analysis to plot results // do not run this code directly

## Reproduction of paper plots:

Paper shows a standard sensitivity analysis over 10,000 runs. 

## Updates to code:

This version includes an update to the vegetation mass calculation to use low latitude rather than global temperautre.

This version also uses a lower 'pre-plant' weathering enhancement than assumed in the version that is published.

This version experiments with adding a direct P source from the oxidaiton of marine DOM during the Shuram excursion.


## Notes:
Requires MATLAB. Tested in R2018a on win10 x64. Run time ~5 seconds. No installation required.


## Modifications for this study (2026)

Modifications © 2026 Ning Lu, Xiamen University. Derived from COPSE V2.1
(Bergman et al., 2004; Lenton et al., 2018; Tostevin & Mills, 2020) and
distributed under the Apache License, Version 2.0 (see `LICENSE`).

Modified files:

- `COPSE_frontend.m` — parameter blocks, initial conditions and switches for
  the marine silica cycle, the alkalinity–DIC closure and the calibration runs.
- `COPSE_equations.m` — state vector expanded from 24 to 27 by adding DSi (25),
  alkalinity (26) and DIC (27); new flux terms for reverse weathering,
  biogenic and abiotic silica burial, carbonate burial and air–sea exchange.

Main changes:

- A marine silica cycle with an explicit dissolved-silica reservoir.
- A temperature-dependent reverse-weathering flux,
  Si_rw = k_rw (DSi/DSi0) exp[Q_rw (T - T_ref)], with Q_rw = 0.0608 K^-1,
  T_ref = 289 K and a flux cap of 15 Tmol Si yr^-1.
- Alkalinity and DIC mass balance coupled to carbonate chemistry, with a
  saturation-dependent carbonate burial flux and air–sea CO2 exchange.

Components adapted from earlier models are credited in the manuscript:
the silica cycle follows Krissansen-Totton & Catling (2020); the alkalinity
layer follows Chen & Lenton (2026); the alkalinity–DIC bookkeeping and
air–sea exchange follow the EaSy/pCO2SYS conventions (Zeebe &
Wolf-Gladrow, 2001; Dickson & Goyet, 1994).

For a file-level diff against the unmodified V2.1 release, compare the
archived version of this folder with the upstream COPSE V2.1 package.

## Sensitivity and figure scripts

The scripts in `sensitivity/` require `chassis` to be set to the path of the
release root and write their results to `outputs/` (created automatically).
Recommended order: main
statistics first (`warmtail_LHS_calibrated.m`, then the two MCMC branches
`warmtail_MCMC.m` and `calibrate_MCMC.m`, then
`fig2_prior_vs_posterior.m`), followed by the remaining figure scripts and,
last, the single-parameter scans.

- `warmtail_LHS_calibrated.m` — random sampling in the observationally admissible region (k_rw × Alk:Si ≤ 2.8) with paired Q_rw = 0 runs; produces the main W_extra distribution.
- `warmtail_MCMC.m` — Metropolis-Hastings posterior of {k_rw, Q_rw, Alk:Si} with a 2-term likelihood (modern pCO2 + DSi) and posterior W_extra.
- `calibrate_MCMC.m` — the same MCMC with a 10-term likelihood (six modern-calibration terms plus four silicon-budget terms).
- `fig2_modern_validation.m` — Fig. 2:modern-state calibration check (model/observation ratios and silicon fluxes versus literature).
- `fig3_paired_runs_figure.m` — composite figure of GMST, DSi/Si_rw, f_T,rw and the geological time bar for paired Q_rw on/off runs.
- `fig4_prior_vs_posterior.m` —Fig. 4: prior propagation versus the 2-term and 10-term posteriors; reads the three `.mat` caches, so run the three scripts above first. 
- `fig5_rw_off.m` — Fig. 5: three-scenario time series (calibrated, no reverse weathering, no temperature dependence).
- `fig6_five_param_wspan.m` — merges five one-at-a-time scans and compares the resulting W_extra ranges.
- `Qrw_sensitivity_scan.m` — scans Q_rw (0–0.12 K^-1), maps it to Ea and Q10, and reports Early Triassic amplification.
- `krw_si_scan.m` — scans k_rw_si over 1.0/1.5/2.0 Tmol Si/yr (estuarine value to Rahman upper bound).
- `AlkSi_ratio_sensitivity.m` — scans the Alk:Si stoichiometry from 1.33 to 4.0 and its effect on W_extra.
- `abiotic_thresh_sensitivity.m` — scans the Isson (2022) abiotic-silica threshold ∂ over 0.6–0.9 mM.
- `Tref_scan_current_chassis.m` — scans the reference temperature T_ref (287–291 K) on the final chassis.
- `sirwmax_sensitivity.m` — scans the reverse-weathering cap Si_rw_max (15–40 Tmol/yr).
- `f_si_evo_PTB_sensitivity.m` — scans the PTB f_si_evo valley depth and reports peak timing (argmax).
- `lowdegassing_pCO2_sensitivity.m` — low-degassing scenario; bisects the degassing scale to reach modern pCO2 = 280 ppm.
- `Wextra_fulltime.m` — full-run ΔT(t) and cumulative W(t) curves for Q_rw on/off.
- `Qrw_sign_symmetry.m` — sign-symmetry control: W_extra for negative and positive Q_rw of equal magnitude.
- `Ksp_taudic_scan.m` — scans K_sp_cal (optionally also k_mccb_omega) against modern pCO2/T/ALK/Ω and the DIC residence time.
- `kcarbw_ksfw_scan.m` — exploratory scan of the base carbon parameters k_carbw/k_sfw and the DIC residence time.
- `warmtail_LHS.m` — Latin hypercube propagation over k_rw × Q_rw × Alk:Si with rejection sampling.
- `substrate_limit_diag.m` — diagnoses substrate self-limitation of reverse weathering, including a cap-removal control.

Note: `Ksp_taudic_scan.m` and `kcarbw_ksfw_scan.m` require the
frontend override hooks documented in their file headers.
