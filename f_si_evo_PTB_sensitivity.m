% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% f_si_evo PTB valley sensitivity scan (with argmax peak timing)
%%%% Science question: a deeper siliceous-plankton collapse (chert gap) → higher Early Triassic DSi →
%%%%           stronger reverse weathering → warmer climate? (tests the collapse→DSi→reverse-weathering→warmth chain)
%%%% Method: fix Q_rw = 0.0608 and scan the f_si_evo valley [0.35 no collapse, 0.10, 0.04, 0.02]
%%%%
%%%% Outputs (two windows):
%%%%   ET window −253 → −248 Ma: mean + peak and its timing (argmax)
%%%%   WIN window −248.7 → −237 Ma (W_extra convention): mean + peak and its timing
%%%%   * the argmax timing tests whether the collapse moves the window maximum from the late part
%%%%     toward the extinction

clear
close all
clc
clear global SENS Gtune

global SENS
global Gtune
Gtune = [];

set(0,'DefaultFigureVisible','off')

% ---- paths (the script can be run from any directory) ----
chassis = '';   % <-- set this to the path of the COPSE-main-RARE root folder
outdir  = fullfile(chassis,'outputs');
if ~exist(outdir,'dir'), mkdir(outdir); end
if isempty(chassis)
    error('Set chassis to the path of the COPSE-main-RARE root folder before running.');
end
cd(chassis)

Q_rw_fixed = 0.0608;                      % fixed temperature dependence (default setting)
valley_list = [0.35, 0.10, 0.04, 0.02];   % valley: 0.35 = no-collapse control; smaller = deeper collapse

et_lo = -253;   et_hi = -248;             % ET window
w0    = -248.7; w1    = -237;             % W_extra clean window

nV = numel(valley_list);

% ET window
DSi_ET   = zeros(nV,1);  Si_rw_ET = zeros(nV,1);
T_ET     = zeros(nV,1);  Tpk_ET   = zeros(nV,1);  tpk_ET = zeros(nV,1);
pCO2_ET  = zeros(nV,1);  pCpk_ET  = zeros(nV,1);

% WIN window
DSi_W    = zeros(nV,1);  Si_rw_W  = zeros(nV,1);
T_W      = zeros(nV,1);  Tpk_W    = zeros(nV,1);  tpk_W  = zeros(nV,1);
pCO2_W   = zeros(nV,1);  pCpk_W   = zeros(nV,1);

for k = 1:nV
    V = valley_list(k);
    fprintf('Running f_si_evo valley = %.2f ...\n', V)

    SENS = struct();
    SENS.Q_rw         = Q_rw_fixed;
    SENS.f_si_evo_PTB = V;
    run = COPSE_frontend(0);
    close all
    s = run.state;

    t  = s.time_myr(:);
    Tc = s.tempC(:);
    Dc = s.DSi_conc_uM(:);
    Sr = s.Si_rw(:)/1e12;
    Pc = s.RCO2(:)*280;

    % ---- ET window ----
    iE = (t >= et_lo) & (t <= et_hi);
    tE = t(iE);  TE = Tc(iE);
    DSi_ET(k)  = mean(Dc(iE));
    Si_rw_ET(k)= mean(Sr(iE));
    T_ET(k)    = mean(TE);
    pCO2_ET(k) = mean(Pc(iE));
    [Tpk_ET(k), jE] = max(TE);
    tpk_ET(k)  = tE(jE);
    pCpk_ET(k) = max(Pc(iE));

    % ---- WIN window ----
    iW = (t >= w0) & (t <= w1);
    tW = t(iW);  TW = Tc(iW);
    DSi_W(k)   = mean(Dc(iW));
    Si_rw_W(k) = mean(Sr(iW));
    T_W(k)     = mean(TW);
    pCO2_W(k)  = mean(Pc(iW));
    [Tpk_W(k), jW] = max(TW);
    tpk_W(k)   = tW(jW);
    pCpk_W(k)  = max(Pc(iW));

    fprintf('VALLEY %.2f | ET Tmean=%.2f Tpk=%.2f@%.1f pCmean=%.0f pCpk=%.0f DSi=%.1f Si_rw=%.2f | WIN Tmean=%.2f Tpk=%.2f@%.1f pCmean=%.0f DSi=%.1f Si_rw=%.2f\n', ...
        V, T_ET(k), Tpk_ET(k), tpk_ET(k), pCO2_ET(k), pCpk_ET(k), DSi_ET(k), Si_rw_ET(k), ...
        T_W(k),  Tpk_W(k),  tpk_W(k),  pCO2_W(k),  DSi_W(k),   Si_rw_W(k));
end

% ---- summary table ----
fprintf('\n===== f_si_evo PTB valley sensitivity (Q_rw = %.4f) =====\n', Q_rw_fixed);
fprintf('%-8s | %-45s | %-45s\n', 'valley', 'ET  (-253 ~ -248 Ma)', 'WIN (-248.7 ~ -237 Ma)');
fprintf('%-8s | %-45s | %-45s\n', '', 'DSi  Si_rw  Tmean  Tpk@t   pCO2', 'DSi  Si_rw  Tmean  Tpk@t   pCO2');
for k = 1:nV
    fprintf('%-8.2f | %5.1f %5.2f %6.2f %6.2f@%-6.1f %5.0f | %5.1f %5.2f %6.2f %6.2f@%-6.1f %5.0f\n', ...
        valley_list(k), DSi_ET(k), Si_rw_ET(k), T_ET(k), Tpk_ET(k), tpk_ET(k), pCO2_ET(k), ...
        DSi_W(k), Si_rw_W(k), T_W(k), Tpk_W(k), tpk_W(k), pCO2_W(k));
end

k35 = find(valley_list==0.35);
fprintf('\n[check] 0.35→0.02: ET ΔTmean = %+.2f °C, WIN ΔTmean = %+.2f °C\n', ...
    T_ET(end)-T_ET(k35), T_W(end)-T_W(k35));
fprintf('[check] 0.35→0.02: WIN peak change = %+.2f °C, peak timing %.1f → %.1f Ma\n', ...
    Tpk_W(end)-Tpk_W(k35), tpk_W(k35), tpk_W(end));
fprintf('[check] 0.35→0.02: ET  peak change = %+.2f °C, peak timing %.1f → %.1f Ma\n', ...
    Tpk_ET(end)-Tpk_ET(k35), tpk_ET(k35), tpk_ET(end));

save(fullfile(outdir,'f_si_evo_PTB_with_argmax.mat'), ...
    'valley_list','DSi_ET','Si_rw_ET','T_ET','Tpk_ET','tpk_ET','pCO2_ET','pCpk_ET', ...
    'DSi_W','Si_rw_W','T_W','Tpk_W','tpk_W','pCO2_W','pCpk_W', ...
    'et_lo','et_hi','w0','w1','Q_rw_fixed');

% ---- figure (four-panel layout retained, peak timings annotated) ----
set(0,'DefaultFigureVisible','on')
figure('Name','f_si_evo PTB valley sensitivity','Color','w','Position',[80 80 1000 700])

subplot(2,2,1)
plot(valley_list, DSi_ET, 'o-', 'LineWidth', 1.8, 'MarkerFaceColor', [0.30 0.45 0.70])
xlabel('f_{si\_evo} valley (smaller = deeper collapse)'); ylabel('Early Triassic DSi (\muM)'); grid on

subplot(2,2,2)
plot(valley_list, Si_rw_ET, 'o-', 'LineWidth', 1.8, 'MarkerFaceColor', [0.47 0.67 0.19])
xlabel('f_{si\_evo} valley'); ylabel('Early Triassic Si_{rw} (Tmol Si/yr)'); grid on

subplot(2,2,3)
plot(valley_list, T_ET, 'o-', 'LineWidth', 1.8, 'MarkerFaceColor', [0.85 0.33 0.10]); hold on
plot(valley_list, Tpk_ET, 's--', 'LineWidth', 1.4, 'MarkerFaceColor', [0.85 0.55 0.10])
xlabel('f_{si\_evo} valley'); ylabel('Early Triassic temperature (C)'); grid on
legend({'ET mean','ET peak'}, 'Location','northwest')

subplot(2,2,4)
plot(valley_list, T_W, 'o-', 'LineWidth', 1.8, 'MarkerFaceColor', [0.85 0.33 0.10]); hold on
plot(valley_list, Tpk_W, 's--', 'LineWidth', 1.4, 'MarkerFaceColor', [0.85 0.55 0.10])
for k = 1:nV
    text(valley_list(k), Tpk_W(k)+0.08, sprintf('%.1f Ma', tpk_W(k)), ...
        'FontSize', 8, 'HorizontalAlignment','center')
end
xlabel('f_{si\_evo} valley'); ylabel('W_{extra} window temperature (C)'); grid on
legend({'WIN mean','WIN peak (argmax timing annotated)'}, 'Location','northwest')

fprintf('\nResults saved to %s\n', fullfile(outdir,'f_si_evo_PTB_with_argmax.mat'));
