% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Self-coupled substrate limitation diagnosis
%%%%
%%%% Purpose: demonstrate that the reverse-weathering flux is subject to self-coupled substrate limitation:
%%%%       Si_rw consumes DSi, and the DSi drawdown in turn suppresses further Si_rw amplification.
%%%%
%%%% Rationale:
%%%%   Si_rw = k_rw × (DSi/DSi0) × exp(Q_rw·ΔT)
%%%%   Q_rw = 0  → Si_rw = k_rw × (DSi/DSi0)  → Si_rw is strictly proportional to DSi (identity)
%%%%   Q_rw > 0  → Si_rw/DSi = (k_rw/DSi0)·exp(Q_rw·ΔT) varies with temperature
%%%%
%%%% Three runs:
%%%%   A) Q_rw = 0,       cap = 15  —— pure substrate-control reference (Si_rw ∝ DSi)
%%%%   B) Q_rw = 0.0608,  cap = 15  —— nominal run (Early Triassic Si_rw sits at the 15 cap)
%%%%   D) Q_rw = 0.0608,  cap = 40  —— relaxed cap, removing cap masking to expose the true self-limitation
%%%%       (40 Tmol is a sensitivity probe with no literature source, used only to separate the cap effect)
%%%%
%%%% Counterfactuals (analytic, no re-integration):
%%%%   *_counter = min[ k_rw·(DSi_hold/DSi0)·f_T, cap ] —— DSi frozen at each window peak
%%%%   gap = *_counter − *_real  —— amplification cut off by DSi self-consumption
%%%%
%%%% Output: 3×3 figure + console diagnostics + substrate_limit_diag.mat

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

% ---- common chassis (nominal converged values) ----
krw      = 1.0e12;
krw_Tmol = krw/1e12;      % = 1.0 Tmol/yr
alkrw    = 2.66;
DSi0_uM  = 50;            % modern reference concentration (μM), corresponding to pars.DSi0
QQ       = 0.0608;        % nominal temperature coefficient
cap15    = 15e12;         % nominal cap
cap40    = 40e12;         % sensitivity-probe cap (only to remove cap masking; not a nominal value)

% ---- run A: Q_rw = 0, cap = 15 ----
SENS = struct();
SENS.k_rw_si = krw; SENS.Q_rw = 0.0; SENS.Alk_Si_RW_ratio = alkrw; SENS.Si_rw_max = cap15;
runA = COPSE_frontend(0); close all
sA = runA.state;

% ---- run B: Q_rw = 0.0608, cap = 15 (nominal) ----
SENS = struct();
SENS.k_rw_si = krw; SENS.Q_rw = QQ; SENS.Alk_Si_RW_ratio = alkrw; SENS.Si_rw_max = cap15;
runB = COPSE_frontend(0); close all
sB = runB.state;

% ---- run D: Q_rw = 0.0608, cap = 40 (cap removed) ----
SENS = struct();
SENS.k_rw_si = krw; SENS.Q_rw = QQ; SENS.Alk_Si_RW_ratio = alkrw; SENS.Si_rw_max = cap40;
runD = COPSE_frontend(0); close all
sD = runD.state;

% ---- extract time series (column vectors; the temperature factor is read from s.f_T_rw) ----
t    = sA.time_myr(:);
DSiA = sA.DSi_conc_uM(:);  SrwA = sA.Si_rw(:)/1e12;  TA = sA.tempC(:);  fT_A = sA.f_T_rw(:);
DSiB = sB.DSi_conc_uM(:);  SrwB = sB.Si_rw(:)/1e12;  TB = sB.tempC(:);  fT_B = sB.f_T_rw(:);
DSiD = sD.DSi_conc_uM(:);  SrwD = sD.Si_rw(:)/1e12;  TD = sD.tempC(:);  fT_D = sD.f_T_rw(:);

% length safeguard (use the common length if the three runs have different grids)
n = min([numel(t) numel(DSiA) numel(SrwA) numel(TA) numel(fT_A) ...
         numel(DSiB) numel(SrwB) numel(TB) numel(fT_B) ...
         numel(DSiD) numel(SrwD) numel(TD) numel(fT_D)]);
t=t(1:n); DSiA=DSiA(1:n); SrwA=SrwA(1:n); TA=TA(1:n); fT_A=fT_A(1:n);
DSiB=DSiB(1:n); SrwB=SrwB(1:n); TB=TB(1:n); fT_B=fT_B(1:n);
DSiD=DSiD(1:n); SrwD=SrwD(1:n); TD=TD(1:n); fT_D=fT_D(1:n);

% ---- console diagnostics ----
fprintf('===== Substrate-limitation diagnostics (A: Q=0/15, B: Q=0.0608/15, D: Q=0.0608/40) =====\n')
fprintf('A) Q_rw=0, cap=15 : Si_rw/DSi = %.5f (theory %.5f = k_rw/DSi0, pure substrate control)\n', ...
    mean(SrwA./max(DSiA,eps)), krw_Tmol/DSi0_uM)
fprintf('B) Q_rw=0.0608    : f_T range %.2f to %.2f (full run, varying with temperature)\n', min(fT_B), max(fT_B))

win = (t >= -253) & (t <= -237);
tt  = t(win);

% window peaks for B
[DSipkB, iB] = max(DSiB(win)); tDSiB = tt(iB);
[SrwpkB, jB] = max(SrwB(win)); tSrwB = tt(jB);
% window peaks for D
[DSipkD, iD] = max(DSiD(win)); tDSiD = tt(iD);
[SrwpkD, jD] = max(SrwD(win)); tSrwD = tt(jD);

fprintf('\nEarly Triassic window (-253 to -237 Ma) peaks:\n')
fprintf('  B(cap15) : DSi peak %.1f μM @ %.2f Ma | Si_rw peak %.2f Tmol @ %.2f Ma (peak offset %.2f Myr)\n', ...
    DSipkB, tDSiB, SrwpkB, tSrwB, tSrwB-tDSiB)
fprintf('  D(cap40) : DSi peak %.1f μM @ %.2f Ma | Si_rw peak %.2f Tmol @ %.2f Ma (peak offset %.2f Myr)\n', ...
    DSipkD, tDSiD, SrwpkD, tSrwD, tSrwD-tDSiD)

% correlation: in-window Si_rw–DSi (B and D)
cB = corr(DSiB(win), SrwB(win));
cD_corr = corr(DSiD(win), SrwD(win));
fprintf('\n  in-window Si_rw–DSi correlation: B(cap15) = %.4f; D(cap40) = %.4f\n', cB, cD_corr)

% ================= counterfactual: DSi frozen at each window peak =================
% —— based on B (cap15) ——
SrwB_pred    = min(krw_Tmol*(DSiB/DSi0_uM)     .* fT_B, cap15/1e12);
SrwB_counter = min(krw_Tmol*(DSipkB/DSi0_uM)   *  fT_B, cap15/1e12);
gapB = SrwB_counter - SrwB;
% —— based on D (cap40) ——
SrwD_pred    = min(krw_Tmol*(DSiD/DSi0_uM)     .* fT_D, cap40/1e12);
SrwD_counter = min(krw_Tmol*(DSipkD/DSi0_uM)   *  fT_D, cap40/1e12);
gapD = SrwD_counter - SrwD;

fprintf('\n===== Counterfactuals (DSi frozen at each window peak) =====\n')
fprintf('Consistency check max|SrwB_pred−SrwB| = %.3f; max|SrwD_pred−SrwD| = %.3f Tmol (should be ≈0)\n', ...
    max(abs(SrwB_pred-SrwB)), max(abs(SrwD_pred-SrwD)))
fprintf('B(cap15) frozen at %.1f μM : gap mean %.2f, peak %.2f Tmol | ratio mean %.2f, peak %.2f\n', ...
    DSipkB, mean(gapB(win)), max(gapB(win)), ...
    mean(SrwB_counter(win)./max(SrwB(win),eps)), max(SrwB_counter(win)./max(SrwB(win),eps)))
fprintf('D(cap40) frozen at %.1f μM : gap mean %.2f, peak %.2f Tmol | ratio mean %.2f, peak %.2f\n', ...
    DSipkD, mean(gapD(win)), max(gapD(win)), ...
    mean(SrwD_counter(win)./max(SrwD(win),eps)), max(SrwD_counter(win)./max(SrwD(win),eps)))

% ---- archive first (values survive even if plotting fails) ----
save(fullfile(outdir,'substrate_limit_diag.mat'), 't','DSiA','SrwA','TA','DSiB','SrwB','TB', ...
    'DSiD','SrwD','TD','fT_A','fT_B','fT_D','SrwB_pred','SrwB_counter','gapB', ...
    'SrwD_pred','SrwD_counter','gapD','krw','DSi0_uM','alkrw','QQ','cap15','cap40')
fprintf('\nResults saved to substrate_limit_diag.mat\n')

% ---- plotting (3×3) ----
try
set(0,'DefaultFigureVisible','on')
figure('Name','Substrate-limitation diagnosis (with cap-removal control)','Color','w','Position',[20 20 1560 940])
xl  = [-260 -200];      % overview (Early–Middle Triassic)
xlz = [-253 -237];      % Early Triassic window

% (1,1) A: Q_rw=0, cap15
subplot(3,3,1); cla
yyaxis left;  plot(t, DSiA, 'LineWidth',1.3); ylabel('DSi (μM)')
yyaxis right; plot(t, SrwA, 'LineWidth',1.3); ylabel('Si_{rw} (Tmol/yr)')
xlim(xl); grid on; title('A) Q_{rw}=0, cap=15: pure substrate ∝'); xlabel('time (Ma)')

% (1,2) B: Q_rw=0.0608, cap15
subplot(3,3,2); cla
yyaxis left;  plot(t, DSiB, 'LineWidth',1.3); ylabel('DSi (μM)')
yyaxis right; plot(t, SrwB, 'LineWidth',1.3); ylabel('Si_{rw} (Tmol/yr)')
xlim(xl); grid on; title('B) Q_{rw}=0.0608, cap=15: Si_{rw} capped'); xlabel('time (Ma)')

% (1,3) D: Q_rw=0.0608, cap40
subplot(3,3,3); cla
yyaxis left;  plot(t, DSiD, 'LineWidth',1.3); ylabel('DSi (μM)')
yyaxis right; plot(t, SrwD, 'LineWidth',1.3); ylabel('Si_{rw} (Tmol/yr)')
xlim(xl); grid on; title('D) Q_{rw}=0.0608, cap=40: cap removed'); xlabel('time (Ma)')

% (2,1) Si_rw–DSi phase portrait
subplot(3,3,4); cla
plot(DSiA, SrwA,'.', 'MarkerSize',5); hold on
plot(DSiB, SrwB,'.', 'MarkerSize',5)
plot(DSiD, SrwD,'.', 'MarkerSize',5)
xlabel('DSi (μM)'); ylabel('Si_{rw} (Tmol/yr)'); grid on
legend('A Q=0/15','B 0.0608/15','D 0.0608/40','Location','northwest')
title('Si_{rw}–DSi phase portrait (D should show a retracing arc)')

% (2,2) Early Triassic window for D: does Si_rw decline after the DSi peak? (key evidence for self-limitation)
subplot(3,3,5); cla
yyaxis left;  plot(t, DSiD, 'LineWidth',1.6); hold on; ylabel('DSi (μM)')
yyaxis right; plot(t, SrwD, 'LineWidth',1.6); ylabel('Si_{rw} (Tmol/yr)')
xline(tDSiD,'--r',sprintf('DSi peak %.1f Ma',tDSiD));
xline(tSrwD,':b',sprintf('Si_rw peak %.1f Ma',tSrwD));
xlim(xlz); grid on; title('D) Early Triassic window: Si_{rw} peaks before DSi → substrate self-limitation'); xlabel('time (Ma)')

% (2,3) Early Triassic window: Si_rw comparison, B vs D
subplot(3,3,6); cla
plot(t, SrwB, 'LineWidth',1.4); hold on
plot(t, SrwD, 'LineWidth',1.4)
yline_seg = @(v) plot(xlz, [v v], ':k', 'LineWidth',1.0);
yline_seg(cap15/1e12);
xlim(xlz); grid on
legend('B cap=15','D cap=40','15 Tmol cap','Location','northwest')
xlabel('time (Ma)'); ylabel('Si_{rw} (Tmol/yr)'); title('Cap masking: B is flattened, D exposes the true peak')

% (3,1) counterfactual gap_B
subplot(3,3,7); cla
yyaxis left;  plot(t, gapB, 'LineWidth',1.5); ylabel('gap_B (Tmol/yr)')
yyaxis right; plot(t, TB, 'LineWidth',1.0); ylabel('T (°C)')
xlim(xlz); grid on; title('gap_B = SrwB_{counter} − SrwB (compressed under the cap)'); xlabel('time (Ma)')

% (3,2) counterfactual gap_D
subplot(3,3,8); cla
yyaxis left;  plot(t, gapD, 'LineWidth',1.5); ylabel('gap_D (Tmol/yr)')
yyaxis right; plot(t, TD, 'LineWidth',1.0); ylabel('T (°C)')
xlim(xlz); grid on; title('gap_D = SrwD_{counter} − SrwD (true self-limitation without the cap)'); xlabel('time (Ma)')

% (3,3) temperature trajectories, B vs D
subplot(3,3,9); cla
plot(t, TB, 'LineWidth',1.3); hold on; plot(t, TD, '--', 'LineWidth',1.2)
xlim(xlz); grid on; xlabel('time (Ma)'); ylabel('T (°C)')
legend('B cap=15','D cap=40','Location','best'); title('Temperature trajectories (the two caps nearly coincide)')

catch ME
    fprintf('\nPlotting failed (%s) —— console diagnostics and the .mat file were saved and are unaffected\n', ME.message)
end