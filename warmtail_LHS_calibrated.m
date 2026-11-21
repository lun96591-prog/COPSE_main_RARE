% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Residual-temperature integral: random sampling within the observationally
%%%% admissible region (n = 2000; values reported in the main text). Purpose: draw **purely random**
%%%% samples inside this region (alkalinity-sink strength k_rw × Alk:Si ≤ 2.8, equivalently
%%%% modern pCO2 ≤ 350 ppm), n = 2000, yielding a robust median and 5%–95% / 2.5%–97.5% quantiles
%%%% of W_extra (the upper tail converges as well).
%%%%
%%%% Differences from warmtail_LHS.m:
%%%%   (i) purely random (rand Monte Carlo) rather than Latin hypercube —— with large n, stratification is unnecessary;
%%%%   (ii) the constraint k_rw × Alk:Si ≤ 2.8 is imposed at the sampling stage (rejection sampling),
%%%%      so all 2000 samples fall inside the observationally admissible region;
%%%%   (iii) the unconstrained prior is not reported (it is diagnostic only; see §9.5).
%%%%
%%%% Parameter priors:
%%%%   k_rw ∈ [0.5, 3.0] Tmol/yr; Q_rw ∈ [0, 0.12]; Alk:Si ∈ [1.33, 4.0]
%%%% Physical constraint: k_rw(Tmol) × Alk:Si ≤ 2.8
%%%% Si_rw_max is fixed at 15 Tmol/yr (nominal value).

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
rng(20261017,'twister')     % reproducible

% ---- parameters ----
n        = 2000;             % number of admissible samples (main-text value; upper tail converged)
                         % cost: ~2.9 s per integration, two runs per sample → 2000 samples ≈ 4000 runs ≈ 3.2 h;
                         %       for a quick preliminary check, run a smaller n first (e.g. 400).
p        = 3;
lb       = [0.5e12, 0.0,   1.33];
ub       = [3.0e12, 0.12,  4.00];
Srw_max  = 15e12;
t_start  = -248.7;
t_end    = -237;
maxSink  = 2.8;              % upper bound on k_rw(Tmol) × Alk:Si (observational constraint)

% ---- purely random rejection sampling (points drawn only inside the admissible region) ----
Params = zeros(n, p);
got = 0;
while got < n
    Pb = lb + rand(1000, p) .* (ub - lb);        % 1000 candidates per pass
    keep = (Pb(:,1)/1e12 .* Pb(:,3)) <= maxSink; % alkalinity-sink strength constraint
    new = Pb(keep,:);
    take = min(n - got, size(new,1));
    Params(got+1 : got+take, :) = new(1:take,:);
    got = got + take;
end
fprintf('Samples collected inside the admissible region (k_rw × Alk:Si ≤ %.1f): %d\n', maxSink, n)

% ---- main loop ----
W      = nan(1,n);
pCO2m  = nan(1,n);
ok     = false(1,n);

for k = 1:n
    krw   = Params(k,1);
    Q     = Params(k,2);
    alkrw = Params(k,3);
    if mod(k,50)==0, fprintf('  %d/%d\n', k, n); end
    try
        % sample
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
        pCO2m(k) = mean(run.state.RCO2(idxm))*280;

        % paired baseline (Q_rw = 0)
        SENS.Q_rw = 0.0;
        runp = COPSE_frontend(0);
        close all
        tp = runp.state.time_myr(:);
        Tp = runp.state.tempC(:);
        Tp_int = interp1(tp, Tp, ts, 'linear', 'extrap');

        sel = (ts >= t_start) & (ts <= t_end);
        W(k) = trapz(ts(sel), Tsm(sel) - Tp_int(sel));
        ok(k) = true;
    catch ME
        fprintf('  FAIL k=%d: %s\n', k, ME.message)
    end
end

% ---- statistics ----
Wv = W(ok);
Wsort = sort(Wv);  nv = numel(Wsort);
Wlo  = Wsort(max(1, round(0.025*nv)));   % 2.5%
Wlo5 = Wsort(max(1, round(0.05*nv)));    % 5%
Wmed = median(Wsort);
Whi95= Wsort(min(nv, round(0.95*nv)));   % 95%
Whi  = Wsort(min(nv, round(0.975*nv)));  % 97.5%
Wmean= mean(Wsort);

fprintf('\n===== Residual-temperature integral (admissible region, n = %d) =====\n', n)
fprintf('Succeeded %d, failed %d\n', sum(ok), n-sum(ok))
fprintf('mean %.1f | median %.1f | 5%%-95%% = [%.1f, %.1f] | 2.5%%-97.5%% = [%.1f, %.1f] °C·Myr\n', ...
    Wmean, Wmed, Wlo5, Whi95, Wlo, Whi)
fprintf('Modern self-check (nearly all admissible runs should be ≤350): pCO2 ∈ [%.0f, %.0f] ppm, count >350 = %d\n', ...
    min(pCO2m(ok)), max(pCO2m(ok)), sum(pCO2m(ok)>350))
fprintf('All runs positive: %d\n', all(Wv>0))

% ---- plotting (all figure text in English, ready for the paper/Supplementary) ----
set(0,'DefaultFigureVisible','on')
figure('Name','W_extra — observationally admissible region','Color','w','Position',[60 60 900 380])
histogram(Wv, 30, 'FaceColor',[0.35 0.55 0.85]); hold on
xline(Wmed,'--r',sprintf('median = %.1f',Wmed),'LineWidth',1.5);
xline(Wlo5,'--k','5th pct'); xline(Whi95,'--k','95th pct','LineWidth',1.2)
xlabel('W_{extra} (°C·Myr)'); ylabel('Sample count'); grid on
title(sprintf('W_{extra} in the observationally admissible region (n = %d; 5%%–95%% = [%.1f, %.1f] °C·Myr)', n, Wlo5, Whi95))

save(fullfile(outdir,'warmtail_LHS_calibrated.mat'), 'Params','W','ok','pCO2m', ...
    'Wmean','Wmed','Wlo','Wlo5','Whi95','Whi','t_start','t_end','lb','ub','maxSink')
fprintf('\nResults saved to warmtail_LHS_calibrated.mat\n')