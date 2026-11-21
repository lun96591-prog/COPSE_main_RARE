% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% K_sp_cal scan: matching the four modern terms + DIC residence time
%%%% Purpose: find K_sp_cal (with optional k_mccb_omega) so that modern pCO2/T/ALK/Ω are all matched,
%%%%       and the DIC residence time τ_DIC falls in the target range (toward ~10^5 yr)
%%%% Depends on COPSE_frontend / COPSE_equations from this release.
%%%%
%%%% =============== frontend patch required before running (otherwise K_sp_cal has no effect) ===============
%%%% After the line "pars.k_mccb_omega = ..." in COPSE_frontend.m, add:
%%%%     if ~isempty(SENS)
%%%%         if isfield(SENS,'K_sp_cal'),     pars.K_sp_cal     = SENS.K_sp_cal;     end
%%%%         if isfield(SENS,'k_mccb_omega'), pars.k_mccb_omega = SENS.k_mccb_omega; end
%%%%     end
%%%% Reason: the present SENS block sits in the silicon-cycle section (K_sp_cal undefined) and is overwritten
%%%% by the alkalinity-section defaults;
%%%%       k_mccb_omega is not in the override list at all.
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

Q_rw = 0.0608;                       % fixed reverse-weathering temperature coefficient ( baseline)
Ksp_list = [6e-7, 7e-7, 8.4e-7, 1.0e-6, 1.15e-6, 1.3e-6];

% optional: also scan k_mccb_omega (requires the frontend override hook first)
% usage: set mccb_scale below to a vector to enter a 2D scan
%   mccb_scale = [0.7, 0.85, 1.0, 1.15, 1.3]  (multiples of the default 21.25e12)
mccb_scale = 1.0;                    % default scans K_sp_cal only (1D); a vector gives 2D

win_lo = -1;   win_hi = 0;           % modern window (last 1 Myr)

% ---- target ranges (hit criteria) ----
tgt_pCO2 = [280 430];   % ppm
tgt_ALK  = [2.0 2.6];   % mmol/kg
tgt_Omg  = [1.5 3.0];   % dimensionless
tgt_tau  = [70 130];    % kyr (toward the classic ~10^5 yr)

Nk = length(Ksp_list);
Nm = length(mccb_scale);
ntot = Nk * Nm;

% result columns: Ksp, mccb×, pCO2(ppm), T(C), ALK(mM), Ω, DIC(mM), τ_DIC(kyr), hit count
R = zeros(ntot, 9);

row = 0;
for i = 1:Nk
  for j = 1:Nm
    row = row + 1;
    Ksp = Ksp_list(i);
    SENS = struct();
    SENS.Q_rw = Q_rw;
    SENS.K_sp_cal = Ksp;                       % effective only after the frontend patch
    if Nm > 1
        SENS.k_mccb_omega = 21.25e12 * mccb_scale(j);  % effective only after the frontend patch
    end

    run = COPSE_frontend(0);
    close all

    % ---- self-check: does the override really take effect? ----
    if abs(run.pars.K_sp_cal - Ksp)/Ksp > 1e-6
        warning('[self-check failed] K_sp_cal override not applied; actual = %.2e; patch frontend as described', run.pars.K_sp_cal)
    end

    s = run.state;
    idx = (s.time_myr >= win_lo) & (s.time_myr <= win_hi);

    pCO2 = mean(s.RCO2(idx)) * 280;            % atmospheric pCO2 (ppm)
    T    = mean(s.tempC(idx));                 % °C
    ALK  = mean(s.Alk_conc_mM(idx));           % mmol/kg
    Om   = mean(s.Omega_cal(idx));             % Ω_cal
    DICm = mean(s.DIC_conc_mM(idx));           % mmol/kg

    DICtot = mean(s.DIC(idx));                 % mol C
    input  = mean(2*s.carbw(idx) + 2*s.silw(idx) + s.DOC_ox(idx));  % mol C/yr
    tau_DIC = DICtot / input / 1e3;            % kyr

    hit = [ pCO2>=tgt_pCO2(1) && pCO2<=tgt_pCO2(2), ...
            ALK>=tgt_ALK(1)   && ALK <=tgt_ALK(2), ...
            Om >=tgt_Omg(1)   && Om <=tgt_Omg(2), ...
            tau_DIC>=tgt_tau(1) && tau_DIC<=tgt_tau(2) ];
    nhit = sum(hit);

    if Nm == 1
        mcci = 1.0;
    else
        mcci = mccb_scale(j);
    end

    R(row,:) = [Ksp, mcci, pCO2, T, ALK, Om, DICm, tau_DIC, nhit];

    fprintf('  [%d/%d] K_sp=%.2e done\n', row, ntot, Ksp)
  end
end

% ---- print the full table once ----
fprintf('\nK_sp_cal scan (Q_rw = %.4f)\n', Q_rw)
fprintf('%s\n', repmat('-',1,72))
fprintf('%-9s %-9s %7s %5s %6s %5s %5s %9s | hits\n', ...
    'K_sp_cal','mccb×','pCO2','T(C)','ALK','Ω','DIC','τ_DIC(kyr)')
fprintf('%-9s %-9s %7s %5s %6s %5s %5s %9s | ----\n', ...
    '','','(ppm)','','(mM)','','(mM)','')
for r = 1:ntot
    fprintf('%-9.2e %-9.2f %7.0f %5.2f %6.2f %5.2f %5.2f %9.1f | %d/4\n', ...
        R(r,1),R(r,2),R(r,3),R(r,4),R(r,5),R(r,6),R(r,7),R(r,8),R(r,9))
end
fprintf('%s\n', repmat('-',1,72))
fprintf('Targets: pCO2∈[%g,%g]ppm  ALK∈[%g,%g]mM  Ω∈[%g,%g]  τ_DIC∈[%g,%g]kyr\n', ...
    tgt_pCO2(1),tgt_pCO2(2), tgt_ALK(1),tgt_ALK(2), tgt_Omg(1),tgt_Omg(2), tgt_tau(1),tgt_tau(2))
fprintf('Reference:  converged chassis K_sp_cal = 8.4e-7 → pCO2 ≈327 ppm, T ≈16.1 °C, ALK ≈2.17 mM, Ω ≈2.0, τ_DIC ≈58 kyr\n')
fprintf('Interpretation: K_sp_cal ↑ → Ω↓ → mccb↓ → more DIC held in the ocean → τ_DIC↑ and ALK↑;\n')
fprintf('      but watch whether pCO2 drifts out of the target range (the three are coupled; no free lunch).\n')