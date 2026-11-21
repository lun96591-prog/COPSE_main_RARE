% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% LHS credible interval for the residual-temperature integral
%%%% Purpose: Latin hypercube sampling of three parameters (k_rw × Q_rw × Alk:Si), a paired baseline
%%%%       and residual-temperature integration, yielding the 5%–95% credible interval of W_extra.
%%%%
%%%% Key: use the clean window (Early + Middle Triassic, -248.7 → -237 Ma) to avoid the Late Triassic
%%%%       degassing artefact (96% of the -100 full-window value is an artefact and was discarded;
%%%%       see §9.5 of the results summary).
%%%%
%%%% Design: two runs per sample (the sample plus a paired run with identical parameters and Q_rw = 0);
%%%%       W_extra = ∫(T_sample − T_paired) dt, isolating the pure temperature-dependence contribution.
%%%%
0%%%% Numerical notes:
%%%%   - each sample = two full deep-time integrations; n = 2000 needs ~4000 integrations at ≈2.9 s each → ~3.2 h.
%%%%   - inlined LHS (no Statistics Toolbox dependency) + clear global to avoid contamination + interp1 alignment.

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
rng(42,'twister')

% ---- LHS settings ----
n  = 2000;                       % sample count (matching warmtail_LHS_calibrated, the main admissible-region result)
p  = 3;
lb = [0.5e12, 0.0,   1.33];     % k_rw, Q_rw, Alk:Si
ub = [3.0e12, 0.12,  4.00];

t_start = -248.7;               % integration start (peak)
t_end   = -237;                 % clean-window end (Early + Middle Triassic, avoiding the Late Triassic peak)

Srw_max = 15e12;                % uniform nominal value of 15 Tmol (numerically stable cap)

% ---- inlined LHS + alkalinity-sink strength constraint (rejection sampling) ----
% prior physical constraint: modern pCO2 ≤350 ppm ⟺ alkalinity-sink strength k_rw × Alk:Si ≲ 2.8 Teq/yr
% rejected at the sampling stage → no wasted deep-time integrations; final samples are uniform in the admissible region
maxSink = 2.8;                      % upper bound on k_rw(Tmol) × Alk:Si (set by the modern pCO2 calibration)
Params = zeros(0, p);
while size(Params,1) < n
    m  = max(n, 50);
    Xb = zeros(m, p);
    for j = 1:p
        edges = (0:m)/m;
        mid   = (edges(1:end-1) + edges(2:end))/2;
        jit   = (rand(1,m) - 0.5)/m;
        Xb(:,j) = (mid(randperm(m)) + jit).';
        Xb(:,j) = max(0, min(1, Xb(:,j)));
    end
    Pb   = lb + Xb .* (ub - lb);
    keep = (Pb(:,1)/1e12 .* Pb(:,3)) <= maxSink;    % alkalinity-sink strength constraint
    Params = [Params; Pb(keep,:)];
end
Params = Params(1:n,:);
fprintf('Sampling: alkalinity-sink constraint k_rw × Alk:Si ≤ %.1f —— %d admissible samples drawn\n', maxSink, n)

% ---- main loop ----
W  = nan(1,n);
ok = false(1,n);
pCO2_mod = nan(1,n);    % modern pCO2 for the calibration-consistency filter (≤350 ppm)

for k = 1:n
    krw   = Params(k,1);
    Q     = Params(k,2);
    alkrw = Params(k,3);
    fprintf('[%2d/%d] k_rw=%.2f Q_rw=%.4f Alk:Si=%.2f ... ', ...
        k, n, krw/1e12, Q, alkrw)

    try
        % (1) sample
        SENS = struct();
        SENS.k_rw_si        = krw;
        SENS.Q_rw           = Q;
        SENS.Alk_Si_RW_ratio = alkrw;
        SENS.Si_rw_max      = Srw_max;
        run = COPSE_frontend(0);
        close all
        ts  = run.state.time_myr(:);
        Tsm = run.state.tempC(:);
        idxm = (ts >= -1) & (ts <= 0);
        pCO2_mod(k) = mean(run.state.RCO2(idxm))*280;   % modern pCO2 (ppm)

        % (2) paired baseline (same k_rw/Alk:Si, Q_rw = 0)
        SENS = struct();
        SENS.k_rw_si        = krw;
        SENS.Q_rw           = 0.0;
        SENS.Alk_Si_RW_ratio = alkrw;
        SENS.Si_rw_max      = Srw_max;
        runp = COPSE_frontend(0);
        close all
        tp = runp.state.time_myr(:);
        Tp = runp.state.tempC(:);

        % align the time grids
        Tp_int = interp1(tp, Tp, ts, 'linear', 'extrap');

        % clean-window residual-temperature integral
        sel = (ts >= t_start) & (ts <= t_end);
        W(k) = trapz(ts(sel), Tsm(sel) - Tp_int(sel));

        ok(k) = true;
        fprintf('W=%.1f\n', W(k))
    catch ME
        fprintf('FAIL: %s\n', ME.message)
    end
end

% ---- quantile statistics (hand-written, no toolbox) ----
% (1) full prior set (no calibration filter)
Wv = W(ok);
Wsort = sort(Wv);  nv = numel(Wsort);
Wlo  = Wsort(max(1, round(0.05*nv)));
Wmed = median(Wsort);
Whi  = Wsort(min(nv, round(0.95*nv)));

% (2) calibration-consistent subset: keep only samples with modern pCO2 ≤ 350 ppm
%    (the modern silicon budget and pCO2 are hard observational constraints; inconsistent parameter
%    sets must not enter the conclusions)
pass = ok & (pCO2_mod <= 350);
Wp = W(pass);  Wp_sort = sort(Wp);  np = numel(Wp_sort);
Wlo_p  = Wp_sort(max(1, round(0.05*np)));
Wmed_p = median(Wp_sort);
Whi_p  = Wp_sort(min(np, round(0.95*np)));

fprintf('\n===== LHS residual-temperature integral (clean window -248.7→-237, cap = 15) =====\n')
fprintf('Samples n = %d, succeeded %d, rejected %d\n', n, sum(ok), n-sum(ok))
fprintf('(1) admissible sampling region (k_rw × Alk:Si ≤ 2.8 imposed): median %.1f | 5%%-95%% = [%.1f, %.1f] °C·Myr\n', ...
    Wmed, Wlo, Whi)
if np > 5
    fprintf('(2) calibration-consistent subset (modern pCO2 ≤350): %d samples (%.0f%%) | median %.1f | 5%%-95%% = [%.1f, %.1f] °C·Myr\n', ...
        np, 100*np/max(1,sum(ok)), Wmed_p, Wlo_p, Whi_p)
    fprintf('   All positive: prior %d | subset %d\n', all(Wv>0), all(Wp>0))
else
    fprintf('(2) the calibration-consistent subset is too small (%d) —— check whether the prior is too wide\n', np)
end

% ---- plotting ----
set(0,'DefaultFigureVisible','on')
figure('Name','Residual-temperature integral LHS','Color','w','Position',[50 50 1000 400])

subplot(1,2,1)
histogram(Wv, 15); hold on
xline(Wmed,'--r',sprintf('median=%.1f',Wmed));
xline(Wlo,'--k'); xline(Whi,'--k')
xlabel('W_{extra} (°C·Myr)'); ylabel('Sample count'); grid on
title(sprintf('(1) full prior set  [%.1f, %.1f]', Wlo, Whi))

subplot(1,2,2)
if np > 5
    histogram(Wp, 15); hold on
    xline(Wmed_p,'--r',sprintf('median=%.1f',Wmed_p));
    xline(Wlo_p,'--k'); xline(Whi_p,'--k')
    title(sprintf('(2) calibration-consistent (pCO_2 ≤350)  [%.1f, %.1f]', Wlo_p, Whi_p))
else
    text(0.5,0.5,'subset too small','HorizontalAlignment','center')
end
xlabel('W_{extra} (°C·Myr)'); ylabel('Sample count'); grid on

save(fullfile(outdir,'warmtail_LHS.mat'), 'Params', 'W', 'ok', 'pCO2_mod', 'pass', ...
    'Wmed','Wlo','Whi','Wmed_p','Wlo_p','Whi_p', 't_start', 't_end', 'lb', 'ub')
fprintf('\nResults saved to warmtail_LHS.mat\n')