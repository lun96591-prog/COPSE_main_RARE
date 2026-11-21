% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% MCMC posterior validation of the residual-temperature integral W_extra
%%%%
%%%% Purpose: address the reviewer concern that only prior propagation, not a posterior, was done.
%%%%      Sample {k_rw, Q_rw, Alk:Si} with Metropolis-Hastings,
%%%%      using a likelihood built on modern-calibration observations (modern pCO2 and modern DSi),
%%%%      then compute W_extra for the posterior samples to obtain the posterior median and 5%–95% quantiles.
%%%%
%%%% Key honesty point (must be stated in the paper/Methods):
%%%%   W_extra = ∫(T_Qrw − T_Q0) dt is counterfactual, and no observation directly constrains it,
%%%%   so no "likelihood for W_extra" exists. What can be done is a modern-calibration likelihood:
%%%%   the model must match modern observations to be accepted, and W_extra is reported as a derived posterior quantity.
%%%%   If the posterior interval ≈ the prior (rejection) interval [2.1, 38.6], prior propagation is sufficient.
%%%%
%%%% Numerical cost:
%%%%   each proposal needs at least one deep-time integration (for modern pCO2/DSi);
%%%%   each accepted sample needs one paired run with Q_rw = 0 to compute W_extra.
%%%%   N_iter = 2000 with ~30% acceptance ≈ 2600 integrations (a few to tens of minutes, machine dependent).
%%%%   Note: fully inlined, with no local functions (to avoid MATLAB scoping/name-resolution issues).

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
rng(20261017,'twister')

% ---- prior bounds and hard constraints ----
lb = [0.5e12, 0.0,   1.33];     % [k_rw, Q_rw, Alk:Si]
ub = [3.0e12, 0.12,  4.00];
maxSink = 2.8;                  % hard upper bound on k_rw(Tmol) × Alk:Si (admissible region)
Srw_max = 15e12;                % nominal reverse-weathering cap
t_start = -248.7;               % clean window
t_end   = -237;

% ---- modern-calibration likelihood (data) ----
pCO2_obs = 280;   pCO2_sig = 50;   % ppm (pre-industrial ~280±10; sigma = 50 absorbs the ~42 ppm
                                  % model bias, so 350 lies at ~1.4 sigma, consistent with the "≤350 tolerance")
DSi_obs  = 80;    DSi_sig  = 20;   % μM (whole-ocean volume mean 70–100)
% Note: pCO2_obs must be the true pre-industrial observation (280), not the model value 321.6
% (otherwise the model fits itself and the likelihood is meaningless)

% ---- MCMC settings ----
N_iter = 2400;        % total iterations; 2400 − burn 400 = 2000 posterior samples (matching prior n = 2000)
                      % 1000 is enough for a trial run; cost ≈ N_iter × 4.4 s ≈ 3 h (DO_WEXTRA = 1)
burn   = 400;         % burn-in
step   = [0.15e12, 0.02, 0.30];   % per-dimension Gaussian random-walk step (tuned to ~25–40% acceptance)

% ---- initial point (main chassis) ----
theta  = [1.0e12, 0.0608, 2.66];
krw = theta(1); Q = theta(2); alkrw = theta(3);

% initial sample (Q_rw = Q): one integration
SENS = struct();
SENS.k_rw_si         = krw;
SENS.Q_rw            = Q;
SENS.Alk_Si_RW_ratio = alkrw;
SENS.Si_rw_max       = Srw_max;
run = COPSE_frontend(0); close all
ts0  = run.state.time_myr(:);
Tsm0 = run.state.tempC(:);
idxm = (ts0 >= -1) & (ts0 <= 0);
pCO2m = mean(run.state.RCO2(idxm))*280;
DSim  = mean(run.state.DSi_conc_uM(idxm));
LL_cur = -0.5*((pCO2m - pCO2_obs)/pCO2_sig)^2 - 0.5*((DSim - DSi_obs)/DSi_sig)^2;

% initial paired run (Q_rw = 0) for W_extra: one integration
SENS.Q_rw = 0.0;
runp = COPSE_frontend(0); close all
tp  = runp.state.time_myr(:);
Tp  = runp.state.tempC(:);
Tp_int = interp1(tp, Tp, ts0, 'linear', 'extrap');
sel0 = (ts0 >= t_start) & (ts0 <= t_end);
W_cur = trapz(ts0(sel0), Tsm0(sel0) - Tp_int(sel0));

if ~isfinite(LL_cur)
    error('Invalid likelihood at the initial point; check that the model runs')
end

% ---- storage ----
accept = 0;
Post   = nan(N_iter, 3);      % posterior samples of theta
Post_W = nan(N_iter, 1);      % posterior W_extra
Npost  = 0;

fprintf('===== MCMC (modern-calibration likelihood) =====\n')
for i = 1:N_iter
    % Gaussian random-walk proposal
    prop = theta + step .* randn(1,3);
    % out of bounds / hard-constraint violation → reject (theta unchanged, likelihood and pair run skipped)
    if ~( any(prop < lb) || any(prop > ub) || (prop(1)/1e12*prop(3) > maxSink) )
        krw = prop(1); Q = prop(2); alkrw = prop(3);
        % proposal sample (one integration)
        SENS = struct();
        SENS.k_rw_si         = krw;
        SENS.Q_rw            = Q;
        SENS.Alk_Si_RW_ratio = alkrw;
        SENS.Si_rw_max       = Srw_max;
        run = COPSE_frontend(0); close all
        ts  = run.state.time_myr(:);
        Tsm = run.state.tempC(:);
        idxm = (ts >= -1) & (ts <= 0);
        pCO2m = mean(run.state.RCO2(idxm))*280;
        DSim  = mean(run.state.DSi_conc_uM(idxm));
        LL_prop = -0.5*((pCO2m - pCO2_obs)/pCO2_sig)^2 - 0.5*((DSim - DSi_obs)/DSi_sig)^2;

        % Metropolis-Hastings acceptance
        if log(rand) < (LL_prop - LL_cur)
            theta  = prop;
            LL_cur = LL_prop;
            accept = accept + 1;
            % run the paired simulation (one integration) only after acceptance, for W_extra
            SENS.Q_rw = 0.0;
            runp = COPSE_frontend(0); close all
            tp  = runp.state.time_myr(:);
            Tp  = runp.state.tempC(:);
            Tp_int = interp1(tp, Tp, ts, 'linear', 'extrap');
            sel = (ts >= t_start) & (ts <= t_end);
            W_cur = trapz(ts(sel), Tsm(sel) - Tp_int(sel));
        end
    end

    % record (after burn-in)
    if i > burn
        Npost = Npost + 1;
        Post(Npost,:) = theta;
        Post_W(Npost) = W_cur;
    end
    if mod(i,500)==0, fprintf('  iter %d/%d | accept %.1f%% | W=%.1f\n', i, N_iter, 100*accept/i, W_cur); end
end

Post   = Post(1:Npost,:);
Post_W = Post_W(1:Npost);

% ---- posterior statistics ----
Wsort = sort(Post_W); nv = numel(Wsort);
Wmed  = median(Wsort);
Wlo5  = Wsort(max(1,round(0.05*nv)));
Whi95 = Wsort(min(nv,round(0.95*nv)));
Wmean = mean(Wsort);

fprintf('\n===== MCMC posterior of the residual-temperature integral (modern-calibration likelihood, n = %d) =====\n', Npost)
fprintf('acceptance %.1f%% | mean %.1f | median %.1f | 5%%-95%% = [%.1f, %.1f] °C·Myr\n', ...
    100*accept/N_iter, Wmean, Wmed, Wlo5, Whi95)
fprintf('Sign positive over the whole posterior: %d\n', all(Post_W>0))
fprintf('Reference: prior propagation in the admissible region (n = 2000), median 19.0 [2.1, 38.6]\n')

% ---- posterior mean of theta ----
fprintf('Posterior mean theta: k_rw = %.2f Tmol | Q_rw = %.4f | Alk:Si = %.2f\n', ...
    mean(Post(:,1))/1e12, mean(Post(:,2)), mean(Post(:,3)))

% ---- plotting ----
set(0,'DefaultFigureVisible','on')
figure('Name','W_extra MCMC posterior','Color','w','Position',[60 60 1100 400])
subplot(1,3,1); histogram(Post_W,25,'FaceColor',[0.35 0.55 0.85]); hold on
xline(Wmed,'--r',sprintf('med=%.1f',Wmed)); xline(Wlo5,'--k'); xline(Whi95,'--k')
xlabel('W_{extra} (°C·Myr)'); ylabel('Sample count'); grid on
title(sprintf('Posterior W_{extra} [%.1f, %.1f]', Wlo5, Whi95))
subplot(1,3,2); scatter(Post(:,1)/1e12, Post_W,8,'filled'); xlabel('k_{rw} (Tmol)'); ylabel('W_{extra}'); grid on; title('k_{rw} vs W')
subplot(1,3,3); scatter(Post(:,2), Post_W,8,'filled'); xlabel('Q_{rw}'); ylabel('W_{extra}'); grid on; title('Q_{rw} vs W')

save(fullfile(outdir,'warmtail_MCMC.mat'), 'Post','Post_W','Wmed','Wlo5','Whi95','Wmean', ...
    'accept','N_iter','burn','lb','ub','maxSink','t_start','t_end','pCO2_obs','DSi_obs')
fprintf('\nResults saved to warmtail_MCMC.mat\n')