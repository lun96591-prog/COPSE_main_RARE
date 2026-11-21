% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Sensitivity of the base carbon parameters k_carbw / k_sfw (exploratory, not a recommended change)
%%%% Purpose: quantify how the COPSE-reloaded base carbon-cycle parameters affect the modern DIC residence time τ_DIC,
%%%%       testing whether the "high weathering input flux" is set by k_carbw / k_sfw (= k_ccdeg − k_silw).
%%%% Note: this is a sensitivity exploration only; changing the base carbon parameters is **not recommended**
%%%% (it would modify the original model rather than extend it).
%%%% Depends on COPSE_frontend / COPSE_equations from this release.
%%%%
%%%% =============== frontend patch required before running (otherwise k_carbw/k_sfw have no effect) ===============
%%%% After the line "pars.k_basw = pars.k_silw * basfrac;" in COPSE_frontend.m, add:
%%%%     %%%% --- SENS override of base carbon fluxes (recompute the carbon mass balance) ---
%%%%     global SENS
%%%%     if ~isempty(SENS) && (isfield(SENS,'k_carbw') || isfield(SENS,'k_sfw'))
%%%%         if isfield(SENS,'k_carbw'), pars.k_carbw = SENS.k_carbw; end
%%%%         if isfield(SENS,'k_sfw'),   pars.k_sfw   = SENS.k_sfw;   end
%%%%         pars.k_mccb  = pars.k_carbw + pars.k_ccdeg - pars.k_sfw;
%%%%         pars.k_silw  = pars.k_mccb  - pars.k_carbw;
%%%%         pars.k_granw = pars.k_silw * (1 - basfrac);
%%%%         pars.k_basw  = pars.k_silw * basfrac;
%%%%     end
%%%% Note: only the carbon chain is overridden; second-order knock-on effects on Sr_sfw/Sr_sedb and P weathering
%%%% are not recomputed (negligible for τ_DIC).
%%%% ═══════════════════════════════════════════════════════════════════════════

clear
close all
clc
% ---- paths ----
chassis = '';   % <-- set this to the path of the COPSE-main-RARE root folder
outdir  = fullfile(chassis,'outputs');
if ~exist(outdir,'dir'), mkdir(outdir); end
if isempty(chassis)
    error('Set chassis to the path of the COPSE-main-RARE root folder before running.');
end
cd(chassis)

global SENS
global Gtune
Gtune = [];

set(0,'DefaultFigureVisible','off')

Q_rw = 0.0608;                        % fixed reverse-weathering temperature coefficient

% scan points (relative to the default carbw = 8e12, sfw = 1.75e12)
k_carbw_list = [6e12, 7e12, 8e12, 9e12, 10e12];   % carbonate weathering (mol C/yr)
k_sfw_list   = [1.0e12, 1.75e12, 2.5e12];          % seafloor weathering (mol C/yr)

win_lo = -1;  win_hi = 0;            % modern window (last 1 Myr)

Nc = length(k_carbw_list);
Ns = length(k_sfw_list);
ntot = Nc * Ns;

% result columns: k_carbw(1e12), k_sfw(1e12), silw(1e12), pCO2, T, ALK(mM), Ω, DIC(mM), τ_DIC(kyr)
R = zeros(ntot, 9);

row = 0;
for i = 1:Nc
  for j = 1:Ns
    row = row + 1;
    kbw = k_carbw_list(i);  ksfw = k_sfw_list(j);

    SENS = struct();
    SENS.Q_rw = Q_rw;
    SENS.k_carbw = kbw;
    SENS.k_sfw   = ksfw;

    run = COPSE_frontend(0);
    close all

    % self-check: override applied (silw should equal 15e12 − ksfw)
    if abs(run.pars.k_carbw - kbw)/kbw > 1e-6 || abs(run.pars.k_sfw - ksfw)/ksfw > 1e-6
        warning('[self-check failed] k_carbw/k_sfw override not applied; patch frontend as described')
    end

    s = run.state;
    idx = (s.time_myr >= win_lo) & (s.time_myr <= win_hi);

    pCO2 = mean(s.RCO2(idx)) * 280;
    T    = mean(s.tempC(idx));
    ALK  = mean(s.Alk_conc_mM(idx));
    Om   = mean(s.Omega_cal(idx));
    DICm = mean(s.DIC_conc_mM(idx));

    silw = mean(s.silw(idx)) / 1e12;                       % Tmol/yr
    DICtot = mean(s.DIC(idx));
    input  = mean(2*s.carbw(idx) + 2*s.silw(idx) + s.DOC_ox(idx));
    tau_DIC = DICtot / input / 1e3;                        % kyr

    R(row,:) = [kbw/1e12, ksfw/1e12, silw, pCO2, T, ALK, Om, DICm, tau_DIC];
    fprintf('  [%d/%d] carbw=%.1f sfw=%.2f done\n', row, ntot, kbw/1e12, ksfw/1e12)
  end
end

fprintf('\nbase carbon-parameter sensitivity (Q_rw = %.4f)\n', Q_rw)
fprintf('%s\n', repmat('-',1,82))
fprintf('%-9s %-7s %-8s %7s %5s %6s %5s %5s %9s\n', ...
    'k_carbw','k_sfw','silw','pCO2','T(C)','ALK','Ω','DIC','τ_DIC(kyr)')
fprintf('%-9s %-7s %-8s %7s %5s %6s %5s %5s %9s\n', ...
    '(1e12)','(1e12)','(1e12)','(ppm)','','(mM)','','(mM)','')
for r = 1:ntot
    fprintf('%-9.1f %-7.2f %-8.2f %7.0f %5.2f %6.2f %5.2f %5.2f %9.1f\n', ...
        R(r,1), R(r,2), R(r,3), R(r,4), R(r,5), R(r,6), R(r,7), R(r,8), R(r,9))
end
fprintf('%s\n', repmat('-',1,82))
fprintf('Expected: silw = 15 − k_sfw (independent of k_carbw);\n')
fprintf('      k_sfw↑ → silw↓ → input↓ → τ_DIC↑; k_carbw↓ → input↓ → τ_DIC↑.\n')
fprintf('Reminder: this explores the base carbonate cycle; changing these original parameters to tune τ is not recommended.\n')