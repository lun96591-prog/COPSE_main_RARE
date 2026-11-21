% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%% save_ensemble_trajectories.m
% Generate Zenodo file 1: 2000 parameter sets × full time trajectories + per-draw W_extra
%
% Window = Early Triassic temperature peak (−248.7 Ma) → end of Middle Triassic (−237 Ma)
%       ≈ 11.7 Myr, corresponding to the integration domain of Eq. (5) in the main text
%
% W_extra(k) = ∫_{-248.7}^{-237} [GMST_active(t) − GMST_counterfactual(t)] dt  (°C·Myr)
%   Computed directly with the trapezoidal rule over the window mask; dt units = Myr
%
% Location: COPSE-main-RARE/sensitivity/  (subfolder)
% Before running, place parameter_samples.csv under COPSE-main-RARE/ (parent directory)
%
% FIX (2026-11-20): ode15s uses adaptive step sizes, so the number of output
% time points varies across runs with different parameters. All state vectors
% are now interpolated onto the reference time axis (from run #1) before
% assignment to the pre-allocated arrays.

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
addpath(chassis)
cd(chassis)

global SENS
global Gtune
Gtune = [];

% ---- read in the 2000 parameter sets ----
T = readtable('parameter_samples.csv');
n = height(T);   % = 2000

% ---- first obtain the full time axis (time unit: Ma, dt ≈ 0.01 Myr) ----
SENS = struct('k_rw_si', T.krw_Tmol(1)*1e12, 'Q_rw', T.Qrw_K(1), 'Alk_Si_RW_ratio', T.Alk_Si(1));
set(0,'DefaultFigureVisible','off')
r0 = COPSE_frontend(0);
close all
t = r0.state.time_myr;
nt = length(t);
dt_myr = mean(abs(diff(t)));   % ≈ 0.01 Myr
fprintf('time axis: nt=%d, t=[%.1f, %.1f] Ma, dt=%.4f Myr\n', nt, min(t), max(t), dt_myr);

% ---- window mask (core: the ≈12 Myr window defined in §2.5 of the main text) ----
win_lo  = -248.7;   % Early Triassic temperature peak
win_hi  = -237.0;   % end of Middle Triassic (integration window)
mask    = (t >= win_lo) & (t <= win_hi);
nwin    = sum(mask);
t_win   = t(mask);
dt_win  = mean(abs(diff(t_win)));   % Myr
fprintf('window : t=[%.1f, %.1f] Ma, nwin=%d, dt=%.4f Myr\n', win_lo, win_hi, nwin, dt_win);

% ---- pre-allocate 3D arrays (nt × n) ----
GMST_a = nan(nt, n); GMST_c = nan(nt, n);
pCO2_a = nan(nt, n); pCO2_c = nan(nt, n);
Si_a   = nan(nt, n); Si_c   = nan(nt, n);
fTa    = nan(nt, n); fTc    = nan(nt, n);
DSi_a  = nan(nt, n); DSi_c  = nan(nt, n);
TA_a   = nan(nt, n); TA_c   = nan(nt, n);
pH_a   = nan(nt, n); pH_c   = nan(nt, n);
Om_a   = nan(nt, n); Om_c   = nan(nt, n);

% ---- per-draw statistics (all computed inside the window) ----
W_extra   = nan(1, n);   % °C·Myr, trapezoidal rule
dT_peak   = nan(1, n);   % °C, active peak − counterfactual peak
dT_mean   = nan(1, n);   % °C, active mean − counterfactual mean
pCO2_modern = nan(1, n); % ppm, last time step of the time axis (modern state)
T_modern    = nan(1, n); % °C
ok_flag     = false(1, n);

tic
fprintf('\n=== Starting %d runs (each draw × 2) ===\n', n);

for k = 1:n
    krw = T.krw_Tmol(k) * 1e12;
    Q   = T.Qrw_K(k);
    alk = T.Alk_Si(k);

    % --- (1) sample run (Q_rw = sampled value) ---
    SENS = struct('k_rw_si', krw, 'Q_rw', Q, 'Alk_Si_RW_ratio', alk);
    ok_a = false;
    try
        ra = COPSE_frontend(0); close all
        s = ra.state;
        ts = s.time_myr(:);   % actual time axis of this run
        GMST_a(:,k) = interp1(ts, s.tempC(:),      t, 'linear', 'extrap');
        pCO2_a(:,k) = interp1(ts, s.RCO2(:)*280,   t, 'linear', 'extrap');
        Si_a(:,k)   = interp1(ts, s.Si_rw(:),      t, 'linear', 'extrap');
        fTa(:,k)    = interp1(ts, s.f_T_rw(:),     t, 'linear', 'extrap');
        DSi_a(:,k)  = interp1(ts, s.DSi(:),        t, 'linear', 'extrap');
        TA_a(:,k)   = interp1(ts, s.ALK(:),        t, 'linear', 'extrap');
        pH_a(:,k)   = interp1(ts, s.pH(:),         t, 'linear', 'extrap');
        Om_a(:,k)   = interp1(ts, s.Omega_cal(:),  t, 'linear', 'extrap');
        ok_a = true;
    catch ME
        fprintf('  sample FAILED k=%d: %s\n', k, ME.message)
    end

    % --- (2) counterfactual run (Q_rw = 0.0) ---
    SENS = struct('k_rw_si', krw, 'Q_rw', 0.0, 'Alk_Si_RW_ratio', alk);
    ok_c = false;
    try
        rc = COPSE_frontend(0); close all
        s = rc.state;
        ts = s.time_myr(:);   % actual time axis of this run
        GMST_c(:,k) = interp1(ts, s.tempC(:),      t, 'linear', 'extrap');
        pCO2_c(:,k) = interp1(ts, s.RCO2(:)*280,   t, 'linear', 'extrap');
        Si_c(:,k)   = interp1(ts, s.Si_rw(:),      t, 'linear', 'extrap');
        fTc(:,k)    = interp1(ts, s.f_T_rw(:),     t, 'linear', 'extrap');
        DSi_c(:,k)  = interp1(ts, s.DSi(:),        t, 'linear', 'extrap');
        TA_c(:,k)   = interp1(ts, s.ALK(:),        t, 'linear', 'extrap');
        pH_c(:,k)   = interp1(ts, s.pH(:),         t, 'linear', 'extrap');
        Om_c(:,k)   = interp1(ts, s.Omega_cal(:),  t, 'linear', 'extrap');
        ok_c = true;
    catch ME
        fprintf('  counterfactual FAILED k=%d: %s\n', k, ME.message)
    end

    ok_flag(k) = ok_a && ok_c;

    % --- (3) W_extra / dT_peak / dT_mean inside the window ---
    if ok_flag(k)
        dT_win = GMST_a(mask,k) - GMST_c(mask,k);
        W_extra(k)   = trapz(t_win, dT_win);   % Myr × °C = °C·Myr
        dT_peak(k)   = max(GMST_a(mask,k)) - max(GMST_c(mask,k));
        dT_mean(k)   = mean(GMST_a(mask,k)) - mean(GMST_c(mask,k));
        pCO2_modern(k) = pCO2_a(end,k);        % last time step
        T_modern(k)    = GMST_a(end,k);
    end

    if mod(k,100)==0
        elapsed = toc;
        rate = k / elapsed;
        eta  = (n - k) / rate;
        fprintf('[k=%4d/%d ok=%4d ok%%=%.0f%% elapsed=%5ds rate=%.1f/s ETA=%5ds]\n', ...
            k, n, sum(ok_flag), sum(ok_flag)/k*100, elapsed, rate, eta);
    end
end

% ---- print summary statistics (corresponds to §3 results in the main text) ----
valid = ok_flag & isfinite(W_extra);
fprintf('\n========================================\n');
fprintf('Window: [%.1f, %.1f] Ma (≈%.1f Myr)\n', win_lo, win_hi, win_hi-win_lo);
fprintf('Successful samples = %d / %d (%.1f%%)\n', sum(valid), n, sum(valid)/n*100);
fprintf('Total elapsed      = %.0f s\n', toc);

fprintf('\nW_extra (°C·Myr):\n');
fprintf('  median       = %.2f\n', median(W_extra(valid)));
fprintf('  5%%–95%%      = [%.2f, %.2f]\n', quantile(W_extra(valid),0.05), quantile(W_extra(valid),0.95));
fprintf('  2.5%%–97.5%%  = [%.2f, %.2f]\n', quantile(W_extra(valid),0.025), quantile(W_extra(valid),0.975));
fprintf('  mean         = %.2f\n', mean(W_extra(valid)));
fprintf('  W > 0 ratio  = %.1f%%\n', sum(W_extra(valid)>0)/sum(valid)*100);

fprintf('\ndT_peak  (°C) : median = %.3f\n', median(dT_peak(valid)));
fprintf('dT_mean  (°C) : median = %.3f\n', median(dT_mean(valid)));

% ---- save (v7.3 supports >2GB; includes window mask metadata) ----
matfile = fullfile(outdir, 'ensemble_trajectories_raw.mat');
save(matfile, ...
    't','dt_myr','nt','n', ...
    'win_lo','win_hi','mask','t_win','dt_win','nwin', ...
    'GMST_a','GMST_c','pCO2_a','pCO2_c','Si_a','Si_c', ...
    'fTa','fTc','DSi_a','DSi_c','TA_a','TA_c','pH_a','pH_c','Om_a','Om_c', ...
    'W_extra','dT_peak','dT_mean','pCO2_modern','T_modern', ...
    'ok_flag', ...
    '-v7.3');

% ---- additionally save per-draw CSV (directly feeds Zenodo counterfactual_pairs.csv) ----
csvfile = fullfile(outdir, 'ensemble_counterfactual_pairs.csv');
out = table( ...
    (1:n)', W_extra', dT_peak', dT_mean', pCO2_modern', T_modern', ok_flag', ...
    'VariableNames', {'draw','W_extra_degC_Myr','dT_peak_C','dT_mean_C', ...
                      'pCO2_modern_ppm','GMST_modern_C','ok_flag'});
writetable(out, csvfile);

fprintf('\nSave complete:\n');
fprintf('  %s    (raw 3D trajectories)\n', matfile);
fprintf('  %s (per-draw W_extra)\n', csvfile);
