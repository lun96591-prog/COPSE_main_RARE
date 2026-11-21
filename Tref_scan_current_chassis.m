% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Tref_scan_current_chassis.m
%%%%
%%%% Purpose: redo the T_ref (reference temperature) sensitivity scan on the **final current chassis**,
%%%%       resolving the inconsistency between the older Tref_sensitivity.m and
%%%%       Qrw_scan_current_chassis.m for the same-configuration clean-window W_extra:
%%%%         old Tref script    : default (T_ref = 289) → W_extra = 10.9 (it also overrides T_ref and integrates on the sample grid)
%%%%         Q_rw chassis       : default            → W_extra = 10.69 (integration on the baseline grid t0)
%%%%       This script mirrors the chassis definition and integration convention of the Q_rw chassis,
%%%%       replacing only the swept parameter with T_ref. Target: T_ref = 289 should reproduce
%%%%       W_extra = 10.69 (difference <1% means the conventions agree).
%%%%
%%%% Chassis (identical to Qrw_scan_current_chassis.m; must be overridden explicitly):
%%%%   k_rw = 1.0 Tmol/yr, Alk:Si = 2.66, Si_rw_max = 15 Tmol/yr, Q_rw = 0.0608
%%%% Scan: T_ref ∈ {287, 288, 289, 290, 291} K
%%%%
%%%% Prerequisite (otherwise SENS.T_ref is ignored):
%%%%   the live equations must compute f_T_rw = exp(pars.Q_rw*(TEMP_gast - pars.T_ref)),
%%%%   and the frontend must define a default pars.T_ref (289).
%%%%   * If W_extra is identical across settings, T_ref is not honoured (equations still hard-code 288),
%%%%     and equations L357 must be changed from 288 to pars.T_ref before rerunning; check 5 reports this.
%%%%
%%%% Windows (same as the chassis script):
%%%%   ET window      -253 → -248 Ma
%%%%   clean window   -248.7 → -237 Ma (W_extra convention, 11.7 Myr)
%%%%
%%%% Modern self-check: report both the "t = 0 instantaneous" and the "-1→0 window mean" conventions,
%%%%           to align with the canonical text values (t = 0: 321.4 ppm / 15.99 °C / f_T = 0.9997).

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

% ---- final current chassis (must be overridden explicitly; do not rely on frontend defaults) ----
krw   = 1.0e12;      % Tmol Si/yr
alk   = 2.66;
Srwmx = 15e12;
Q_ref = 0.0608;      % reverse-weathering temperature coefficient (default setting)

% ---- scan points ----
Tref_list = [287, 288, 289, 290, 291];
nT = numel(Tref_list);

% ---- windows ----
et_lo = -253;   et_hi = -248;
w0    = -248.7; w1    = -237;

% ---- preallocation ----
W_cln   = zeros(1,nT);   % clean-window W_extra (°C·Myr)
dTm_cln = zeros(1,nT);   % clean-window mean ΔT (= W/window length)
dT_ET   = zeros(1,nT);   % ET-window mean ΔT
Si_ET   = zeros(1,nT);   % ET-window mean Si_rw (Tmol/yr)
% modern self-check: two conventions
fT_w  = zeros(1,nT);  pC_w  = zeros(1,nT);  T_w  = zeros(1,nT);   % window mean (-1→0)
fT_t0 = zeros(1,nT);  pC_t0 = zeros(1,nT);  T_t0 = zeros(1,nT);   % t = 0 instantaneous

% ---- reference run: Q_rw = 0 ----
%   Note: with Q_rw = 0, f_T_rw = exp(0) = 1 independently of T_ref, so the baseline is run once
%   (identical to Qrw_scan_current_chassis.m: no SENS.T_ref, frontend default used).
fprintf('Reference run Q_rw = 0 (chassis k_rw=%.1f, Alk:Si=%.2f, cap=%.0f) ...\n', ...
    krw/1e12, alk, Srwmx/1e12)
SENS = struct();
SENS.k_rw_si         = krw;
SENS.Q_rw            = 0.0;
SENS.Alk_Si_RW_ratio = alk;
SENS.Si_rw_max       = Srwmx;
r0 = COPSE_frontend(0);
close all
t0 = r0.state.time_myr(:);
T0 = r0.state.tempC(:);

sel_et0  = (t0 >= et_lo) & (t0 <= et_hi);
sel_cln0 = (t0 >= w0)    & (t0 <= w1);

T_ET0 = mean(T0(sel_et0));

% ---- per-setting scan ----
for i = 1:nT
    fprintf('[%d/%d] T_ref = %d K ... ', i, nT, Tref_list(i))

    SENS = struct();
    SENS.k_rw_si         = krw;
    SENS.Q_rw            = Q_ref;
    SENS.Alk_Si_RW_ratio = alk;
    SENS.Si_rw_max       = Srwmx;
    SENS.T_ref           = Tref_list(i);     % ← the only swept parameter
    r1 = COPSE_frontend(0);
    close all

    ts   = r1.state.time_myr(:);
    T1r  = r1.state.tempC(:);
    fT1r = r1.state.f_T_rw(:);

    % * key: state variables such as Si_rw / RCO2 / tempC / f_T_rw must be window-averaged on
    %   their own run grids; only the temperature used for differencing is interpolated to t0.
    sel_et_r = (ts >= et_lo) & (ts <= et_hi);
    idxm_r   = (ts >= -1)    & (ts <= 0);
    [~, i0r] = min(abs(ts - 0));             % index of the sample closest to t = 0

    T_ET_r   = mean(T1r(sel_et_r));
    dT_ET(i) = T_ET_r - T_ET0;
    Si_ET(i) = mean(r1.state.Si_rw(sel_et_r)) / 1e12;

    % clean-window integral: interpolate both runs onto the **common grid t0**, then difference
    T1      = interp1(ts, T1r, t0, 'linear', 'extrap');
    dT      = T1 - T0;
    W_cln(i)   = trapz(t0(sel_cln0), dT(sel_cln0));
    dTm_cln(i) = W_cln(i) / (w1 - w0);

    % modern self-check (two conventions)
    fT_w(i)  = mean(fT1r(idxm_r));   pC_w(i)  = mean(r1.state.RCO2(idxm_r)) * 280;   T_w(i)  = mean(T1r(idxm_r));
    fT_t0(i) = fT1r(i0r);            pC_t0(i) = r1.state.RCO2(i0r) * 280;            T_t0(i) = T1r(i0r);

    fprintf('W_extra=%.2f | ΔT_ET=%+.2f | fT(t0)=%.4f\n', W_cln(i), dT_ET(i), fT_t0(i))
end

% ---- report ----
fprintf('\n===== T_ref sensitivity (current chassis: k_rw=%.1f, Alk:Si=%.2f, cap=%.0f, Q_rw=%.4f) =====\n', ...
    krw/1e12, alk, Srwmx/1e12, Q_ref)
fprintf('%-8s %-11s %-12s %-10s %-11s %-9s\n', 'T_ref','W_extra','ΔT_clean','ΔT_ET','Si_rw_ET','fT(t0)')
fprintf('%s\n', repmat('-', 1, 70))
for i = 1:nT
    fprintf('%-8d %-11.2f %+12.2f %+10.2f %-11.2f %-9.4f\n', ...
        Tref_list(i), W_cln(i), dTm_cln(i), dT_ET(i), Si_ET(i), fT_t0(i))
end

i289 = find(Tref_list==289);
fprintf('\n[check 1] default T_ref = 289 → W_extra = %.2f °C·Myr (should match the Q_rw chassis default 10.69, difference <1%%)\n', ...
    W_cln(i289))
fprintf('[check 2] W_extra range across settings: %.2f to %.2f °C·Myr (text §3.3 reports [9.0, 12.3])\n', ...
    min(W_cln), max(W_cln))
fprintf('[check 3] change relative to the 289 setting: %+.1f%% to %+.1f%%\n', ...
    100*(min(W_cln)-W_cln(i289))/W_cln(i289), 100*(max(W_cln)-W_cln(i289))/W_cln(i289))
fprintf('[check 4] all settings positive: %d (1 = yes)\n', all(W_cln>0))
fprintf('[check 5] T_ref honoured (W_extra varies with T_ref): %d (1 = yes; 0 = equations still hard-code 288)\n', ...
    numel(unique(round(W_cln,3))) > 1)

fprintf('\n---- modern self-check (two conventions, for comparison with the canonical text values) ----\n')
fprintf('%-8s | %-22s | %-22s\n', 'T_ref', 'window mean (-1→0)', 't = 0 instantaneous')
fprintf('%-8s | %-6s %-7s %-6s | %-6s %-7s %-6s\n', '', 'fT','pCO2','T','fT','pCO2','T')
for i = 1:nT
    fprintf('%-8d | %-6.3f %-7.0f %-6.2f | %-6.4f %-7.1f %-6.2f\n', ...
        Tref_list(i), fT_w(i), pC_w(i), T_w(i), fT_t0(i), pC_t0(i), T_t0(i))
end
fprintf('\n  Note: the canonical text values use the "t = 0 instantaneous" convention → 321.4 ppm / 15.99 °C / f_T = 0.9997.\n')
fprintf('      the t = 0 column at T_ref = 289 should match; the window-mean column (337 / 16.33) is **not** the text convention and must not be quoted directly.\n')

% ---- plotting ----
set(0,'DefaultFigureVisible','on')
figure('Name','T_ref sensitivity (current chassis)','Color','w','Position',[60 60 1000 400])

subplot(1,2,1)
plot(Tref_list, W_cln, 'o-', 'LineWidth',1.6, 'MarkerSize',7, 'Color',[0 0.4 0.8]); hold on
yline(0,'-k')
xlabel('T_{ref} (K)'); ylabel('W_{extra} (\circC\cdotMyr)'); grid on
title('Clean-window residual integral vs T_{ref}')

subplot(1,2,2)
plot(Tref_list, fT_t0, 's-', 'LineWidth',1.6, 'MarkerSize',7, 'Color',[0.8 0.4 0]); hold on
yline(1,'--k','f_T=1')
xlabel('T_{ref} (K)'); ylabel('modern f_{T,rw} (t=0)'); grid on
title('Modern f_T self-check vs T_{ref}')

save(fullfile(outdir,'Tref_scan_current_chassis.mat'), 'Tref_list','W_cln','dTm_cln','dT_ET','Si_ET', ...
    'fT_w','pC_w','T_w','fT_t0','pC_t0','T_t0','krw','alk','Srwmx','Q_ref','et_lo','et_hi','w0','w1')
fprintf('\nResults saved to Tref_scan_current_chassis.mat\n')
