% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Abiotic-silica threshold ∂ sensitivity scan (Isson 2022 threshold form)
%%%%
%%%% Purpose: test whether the Isson (2022) abiotic-silica threshold ∂, varied over 0.6–0.9 mM,
%%%%       affects (i) the modern baseline (pCO₂ / DSi) and (ii) the Early Triassic clean-window
%%%%
%%%%       residual-temperature integral W_extra. Form (Isson 2022, Eq. 6):
%%%%   Si_abiotic = 0                        , [Si] ≤ ∂
%%%%   Si_abiotic = C · ( [Si]/∂ − 1 )^1.1   , [Si] > ∂
%%%%   the 0.6–0.9 mol m⁻³ range of ∂ comes from amorphous-silica solubility:
%%%%     Siever (1962); Gunnarsson & Arnórsson (2000); scan range from Isson 2022 supplementary tables
%%%%   This project uses the median ∂ = 0.8 mM by default.
%%%%
%%%% Parameter hook: SENS.Si_abiotic_thresh (same name as in rwoff_abiotic_check.m)
%%%% Other parameters fixed at nominal converged values: k_rw = 1.0, Q_rw = 0.0608, Alk:Si = 2.66,
%%%% Si_rw_max = 15
%%%%
%%%% Expectation: the Early Triassic DSi peak is only ~0.33–0.38 mM, below the lower bound of ∂ (0.6 mM),
%%%%       so Phanerozoic Si_abiotic ≡ 0 and W_extra barely changes with ∂ (∂ is insensitive in the Phanerozoic).
%%%%
%%%% Output: console table + abiotic_thresh_sensitivity.mat + figure

clear
close all
clc
clear global SENS Gtune
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
rng(20261019,'twister')

% ---- fixed chassis (nominal converged values) ----
krw      = 1.0e12;
QQ       = 0.0608;
alkrw    = 2.66;
Srw_max  = 15e12;
t_start  = -248.7;
t_end    = -237;

% ---- scanned thresholds (mM → mol/kg, i.e. ×1e-3) ----
th_list  = [0.6 0.7 0.8 0.9];      % Isson 2022 range; 0.8 is the project default
nth      = numel(th_list);

W        = nan(1,nth);      % clean-window residual-temperature integral (°C·Myr)
pCO2m    = nan(1,nth);      % modern pCO₂ (ppm)
DSIm     = nan(1,nth);      % modern DSi (μM)
DSIetpk  = nan(1,nth);      % Early Triassic window DSi peak (μM)
Sabio_max= nan(1,nth);      % full-period abiotic-silica peak (Tmol/yr)
ok       = false(1,nth);

win = [];   % time window established by the first successful run

for k = 1:nth
    th = th_list(k);
    try
        % ---- run A: Q_rw = QQ (nominal) ----
        SENS = struct();
        SENS.k_rw_si           = krw;
        SENS.Q_rw              = QQ;
        SENS.Alk_Si_RW_ratio   = alkrw;
        SENS.Si_rw_max         = Srw_max;
        SENS.Si_abiotic_thresh = th*1e-3;    % mM → mol/kg
        runA = COPSE_frontend(0); close all
        sA = runA.state;

        ts  = sA.time_myr(:);
        Tsm = sA.tempC(:);
        idxm = (ts >= -1) & (ts <= 0);
        pCO2m(k) = mean(sA.RCO2(idxm))*280;
        DSIm(k)  = mean(sA.DSi_conc_uM(idxm));
        if isfield(sA,'Si_abiotic')
            Sabio_max(k) = max(sA.Si_abiotic(:))/1e12;
        end
        selW = (ts >= -253) & (ts <= -237);
        DSIetpk(k) = max(sA.DSi_conc_uM(selW));
        win = selW;

        % ---- run B: paired baseline Q_rw = 0 ----
        SENS.Q_rw = 0.0;
        runB = COPSE_frontend(0); close all
        sB = runB.state;
        tp = sB.time_myr(:);
        Tp = sB.tempC(:);
        Tp_int = interp1(tp, Tp, ts, 'linear', 'extrap');
        sel = (ts >= t_start) & (ts <= t_end);
        W(k) = trapz(ts(sel), Tsm(sel) - Tp_int(sel));

        ok(k) = true;
        fprintf('  ∂ = %.2f mM : ok  (W=%.2f, pCO2=%.0f, DSi=%.0f, ET-DSi peak=%.0f, Sabio_max=%.3f)\n', ...
            th, W(k), pCO2m(k), DSIm(k), DSIetpk(k), Sabio_max(k))
    catch ME
        fprintf('  ∂ = %.2f mM : FAIL (%s)\n', th, ME.message)
    end
end

% ---- summary table ----
fprintf('\n===== Abiotic-silica threshold ∂ sensitivity (Isson 2022 range 0.6–0.9 mM) =====\n')
fprintf('%-8s %-10s %-10s %-12s %-12s %-10s\n', ...
    '∂(mM)','W_extra','pCO2 mod','DSi mod','ET-DSi peak','S_abio peak')
fprintf('%-8s %-10s %-10s %-12s %-12s %-10s\n', ...
    '','(°C·Myr)','(ppm)','(μM)','(μM)','(Tmol/yr)')
for k = 1:nth
    if ok(k)
        fprintf('%-8.2f %-10.2f %-10.0f %-12.0f %-12.0f %-10.3f\n', ...
            th_list(k), W(k), pCO2m(k), DSIm(k), DSIetpk(k), Sabio_max(k))
    end
end
Wv = W(ok);
if ~isempty(Wv)
    fprintf('\nW_extra range [%.2f, %.2f] °C·Myr, maximum deviation from the 0.8 setting %.1f%%\n', ...
        min(Wv), max(Wv), 100*max(abs(Wv - W(th_list==0.8))/max(abs(W(th_list==0.8)),eps)))
end
if all(Sabio_max(ok) < 1e-3)
    fprintf('Full-period abiotic-silica peak ≈ 0 (threshold never triggered in the Phanerozoic for any ∂) → ∂ is insensitive in the Phanerozoic\n')
end

% ---- archive ----
save(fullfile(outdir,'abiotic_thresh_sensitivity.mat'), 'th_list','W','pCO2m','DSIm', ...
    'DSIetpk','Sabio_max','ok','krw','QQ','alkrw','Srw_max','t_start','t_end')
fprintf('\nResults saved to abiotic_thresh_sensitivity.mat\n')

% ---- plotting (all figure text in English) ----
try
set(0,'DefaultFigureVisible','on')
figure('Name','Abiotic-silica threshold sensitivity','Color','w','Position',[60 60 1000 380])

subplot(1,2,1)
plot(th_list(ok), W(ok), '-o', 'LineWidth',1.5); grid on
xlabel('\partial (mM)'); ylabel('W_{extra} (°C·Myr)')
title('Residual-temperature integral vs threshold')

subplot(1,2,2)
yyaxis left
plot(th_list(ok), DSIetpk(ok), '-s', 'LineWidth',1.5); hold on
yline_lo = min(th_list(ok))*1e3;      % lower bound of ∂ (μM)
plot([th_list(ok(1)) th_list(ok(end))], [yline_lo yline_lo], ':k', 'LineWidth',1.2)
ylabel('Early Triassic DSi peak (\muM)')
yyaxis right
plot(th_list(ok), Sabio_max(ok), '-^', 'LineWidth',1.5); ylabel('max Si_{abiotic} (Tmol yr^{-1})')
xlabel('\partial (mM)'); grid on
title('Early Triassic DSi vs threshold (dotted: lower bound of \partial)')
catch ME
    fprintf('\nPlotting failed (%s) —— values and the .mat file were saved\n', ME.message)
end