% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Modern-calibration posterior MCMC (6-variable calibration + silicon budget,
%%%%
%%%% Purpose: run Metropolis-Hastings over {k_rw, Q_rw, Alk:Si} with a likelihood that no longer
%%%%      uses only pCO2 + DSi, but ten Gaussian terms: six modern-calibration terms plus four
%%%%      modern silicon-budget terms, giving a posterior closer to modern observations; W_extra is
%%%%
%%%% Difference from warmtail_MCMC.m: the likelihood expands from 2 to 10 terms (richer, closer to data).
%%%% The parameter set stays {k_rw, Q_rw, Alk:Si}, giving a direct test of whether richer data change the posterior.
%%%%
%%%% Honesty points (must be stated in the paper):
%%%%   - W_extra is counterfactual with no direct observations, so the likelihood can only use
%%%%     modern calibration and the modern silicon budget;
%%%%   - Si_silw / Si_sfw / Si_bio are controlled mainly by the fixed parameters k_bio_si /
%%%%     k_sfw_si / Si_ALK_ratio,
%%%%     which this MCMC does not tune, so their constraint on k_rw / Alk:Si is indirect
%%%%     (via DSi and climate coupling);
%%%%   - pCO2 / Ω carry a structural model offset (model modern ≈321.6 ppm, Ω ≈2.01 vs
%%%%     observed ~280 ppm and ~4.5),
%%%%     so both terms use wide sigma and act as soft constraints only.
%%%%
%%%% Model output fields (verified against the live version): RCO2, DSi_conc_uM, Alk_conc_mM,
%%%%   Omega_cal, tempC, pH, Si_silw, Si_sfw, Si_bio, Si_rw (fluxes in mol/yr; divide by 1e12 for Tmol/yr).
%%%%
%%%% * Cost (measured ~2.9 s per integration, from 2000 warmtail_MCMC iterations ≈ 2.5 h):
%%%%   total ≈ (N_iter + N_accept) × 2.9 s ≈ N_iter × 4.4 s (with DO_WEXTRA = 1)
%%%%   N_iter = 600 → ~45 min; 1200 → ~1.5 h; 2000 → ~2.5 h; with DO_WEXTRA = 0, ≈N_iter × 2.9 s (halved).
%%%%   Note: fully inlined, no local functions (to avoid MATLAB scoping/function-name resolution issues).

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
rng(20261018,'twister')

% ---- prior bounds and hard constraints ----
lb = [0.5e12, 0.0,   1.33];     % [k_rw, Q_rw, Alk:Si]
ub = [3.0e12, 0.12,  4.00];
maxSink = 2.8;                  % hard upper bound on k_rw(Tmol) × Alk:Si (proxy for modern pCO2 ≤350)
Srw_max = 15e12;                % nominal reverse-weathering cap
t_start = -248.7;               % clean window
t_end   = -237;

% ---- 10-term likelihood targets: {observed value, sigma} (modern; units in comments) ----
% (1) six modern-calibration terms (carbon chemistry + climate)
pCO2_obs = 280;   pCO2_sig = 50;    % ppm (pre-industrial; model ≈321.6 offset → wide sigma)
DSi_obs  = 80;    DSi_sig  = 20;    % μM (whole-ocean volume mean 70–100)
TAlk_obs = 2.30;  TAlk_sig = 0.15;  % mmol/kg (surface total alkalinity ~2.2–2.4)
T_obs    = 15.0;  T_sig    = 1.5;   % °C (modern GMST ~15; model ≈16.3)
pH_obs   = 8.05;  pH_sig   = 0.15;  % surface pH ~8.05 (model ≈8.07)
Om_obs   = 3.0;   Om_sig   = 1.5;   % Ω_calcite (model ≈2.01 biased low; true surface ≈4.5 → wide sigma)
% (2) four modern silicon-budget terms (Tmol Si/yr; Tréguer 2016/2021 + Rahman 2017)
Si_silw_obs = 7.3;  Si_silw_sig = 1.5;   % terrigenous silicate-weathering Si input
Si_sfw_obs  = 1.5;  Si_sfw_sig  = 1.0;   % seafloor-weathering Si input (range 0.8–1.9)
Si_bio_obs  = 6.9;  Si_bio_sig  = 1.5;   % net biogenic-silica burial
Si_rw_obs   = 2.0;  Si_rw_sig   = 1.0;   % reverse weathering / authigenic clay (estuarine lower bound 1.5±0.5)

% ---- MCMC settings ----
N_iter = 2150;        % total iterations; 2150 − burn 150 = 2000 posterior samples (matching prior n = 2000)
                      % trial run 600 ≈45 min; cost ≈ N_iter × 4.4 s ≈ 2.7 h
burn   = 150;        % burn-in
step   = [0.15e12, 0.02, 0.30];   % per-dimension Gaussian random-walk step
DO_WEXTRA = 1;       % 1 = run the paired simulation after acceptance for W_extra (cost ×1.5);
                     % 0 = parameters only (cost halved)

% ---- initial point (main chassis) ----
theta  = [1.0e12, 0.0608, 2.66];
krw = theta(1); Q = theta(2); alkrw = theta(3);

% initial sample: one integration
SENS = struct();
SENS.k_rw_si = krw; SENS.Q_rw = Q; SENS.Alk_Si_RW_ratio = alkrw; SENS.Si_rw_max = Srw_max;
run = COPSE_frontend(0); close all
% one-off field check (confirm the live model provides all output fields before a long run)
need = {'time_myr','tempC','RCO2','DSi_conc_uM','Alk_conc_mM','Omega_cal','pH','Si_silw','Si_sfw','Si_bio','Si_rw'};
have = isfield(run.state, need);
if ~all(have)
    fprintf('Missing fields: '); fprintf('%s ', need{~have}); fprintf('\n')
    error('The live model lacks the fields listed above; check workingstate field names in equations')
end
ts0  = run.state.time_myr(:);
Tsm0 = run.state.tempC(:);
idxm = (ts0 >= -1) & (ts0 <= 0);
pCO2 = mean(run.state.RCO2(idxm))*280;
DSi  = mean(run.state.DSi_conc_uM(idxm));
TAlk = mean(run.state.Alk_conc_mM(idxm));
TT   = mean(run.state.tempC(idxm));
pH   = mean(run.state.pH(idxm));
Om   = mean(run.state.Omega_cal(idxm));
Si_silw = mean(run.state.Si_silw(idxm))/1e12;
Si_sfw  = mean(run.state.Si_sfw(idxm))/1e12;
Si_bio  = mean(run.state.Si_bio(idxm))/1e12;
Si_rw   = mean(run.state.Si_rw(idxm))/1e12;
LL_cur = -0.5*( ((pCO2-pCO2_obs)/pCO2_sig)^2 + ((DSi-DSi_obs)/DSi_sig)^2 + ...
                ((TAlk-TAlk_obs)/TAlk_sig)^2 + ((TT-T_obs)/T_sig)^2 + ...
                ((pH-pH_obs)/pH_sig)^2 + ((Om-Om_obs)/Om_sig)^2 + ...
                ((Si_silw-Si_silw_obs)/Si_silw_sig)^2 + ((Si_sfw-Si_sfw_obs)/Si_sfw_sig)^2 + ...
                ((Si_bio-Si_bio_obs)/Si_bio_sig)^2 + ((Si_rw-Si_rw_obs)/Si_rw_sig)^2 );
% initial paired run for W_extra
if DO_WEXTRA
    SENS.Q_rw = 0.0;
    runp = COPSE_frontend(0); close all
    tp = runp.state.time_myr(:); Tp = runp.state.tempC(:);
    Tp_int = interp1(tp, Tp, ts0, 'linear', 'extrap');
    sel0 = (ts0 >= t_start) & (ts0 <= t_end);
    W_cur = trapz(ts0(sel0), Tsm0(sel0) - Tp_int(sel0));
else
    W_cur = nan;
end
if ~isfinite(LL_cur)
    error('Invalid likelihood at the initial point; check that the model runs')
end

% ---- storage ----
accept = 0;
Post   = nan(N_iter, 3);
Post_W = nan(N_iter, 1);
Npost  = 0;

fprintf('===== MCMC (6 calibration + silicon budget, 10-term likelihood; N_iter = %d, DO_WEXTRA = %d) =====\n', N_iter, DO_WEXTRA)
for i = 1:N_iter
    prop = theta + step .* randn(1,3);
    if ~( any(prop < lb) || any(prop > ub) || (prop(1)/1e12*prop(3) > maxSink) )
        krw = prop(1); Q = prop(2); alkrw = prop(3);
        % proposal sample: one integration
        SENS = struct();
        SENS.k_rw_si = krw; SENS.Q_rw = Q; SENS.Alk_Si_RW_ratio = alkrw; SENS.Si_rw_max = Srw_max;
        run = COPSE_frontend(0); close all
        ts  = run.state.time_myr(:);
        Tsm = run.state.tempC(:);
        idxm = (ts >= -1) & (ts <= 0);
        pCO2 = mean(run.state.RCO2(idxm))*280;
        DSi  = mean(run.state.DSi_conc_uM(idxm));
        TAlk = mean(run.state.Alk_conc_mM(idxm));
        TT   = mean(run.state.tempC(idxm));
        pH   = mean(run.state.pH(idxm));
        Om   = mean(run.state.Omega_cal(idxm));
        Si_silw = mean(run.state.Si_silw(idxm))/1e12;
        Si_sfw  = mean(run.state.Si_sfw(idxm))/1e12;
        Si_bio  = mean(run.state.Si_bio(idxm))/1e12;
        Si_rw   = mean(run.state.Si_rw(idxm))/1e12;
        LL_prop = -0.5*( ((pCO2-pCO2_obs)/pCO2_sig)^2 + ((DSi-DSi_obs)/DSi_sig)^2 + ...
                        ((TAlk-TAlk_obs)/TAlk_sig)^2 + ((TT-T_obs)/T_sig)^2 + ...
                        ((pH-pH_obs)/pH_sig)^2 + ((Om-Om_obs)/Om_sig)^2 + ...
                        ((Si_silw-Si_silw_obs)/Si_silw_sig)^2 + ((Si_sfw-Si_sfw_obs)/Si_sfw_sig)^2 + ...
                        ((Si_bio-Si_bio_obs)/Si_bio_sig)^2 + ((Si_rw-Si_rw_obs)/Si_rw_sig)^2 );

        if log(rand) < (LL_prop - LL_cur)
            theta  = prop;
            LL_cur = LL_prop;
            accept = accept + 1;
            % paired run for W_extra after acceptance (optional)
            if DO_WEXTRA
                SENS.Q_rw = 0.0;
                runp = COPSE_frontend(0); close all
                tp = runp.state.time_myr(:); Tp = runp.state.tempC(:);
                Tp_int = interp1(tp, Tp, ts, 'linear', 'extrap');
                sel = (ts >= t_start) & (ts <= t_end);
                W_cur = trapz(ts(sel), Tsm(sel) - Tp_int(sel));
            else
                W_cur = nan;
            end
        end
    end

    if i > burn
        Npost = Npost + 1;
        Post(Npost,:) = theta;
        Post_W(Npost) = W_cur;
    end
    if mod(i,200)==0, fprintf('  iter %d/%d | accept %.1f%% | W=%.1f\n', i, N_iter, 100*accept/i, W_cur); end
end

Post   = Post(1:Npost,:);
Post_W = Post_W(1:Npost);

% ---- posterior statistics (with DO_WEXTRA = 0, W is all NaN and statistics are skipped) ----
if DO_WEXTRA
    Wsort = sort(Post_W); nv = numel(Wsort);
    Wmed  = median(Wsort);
    Wlo5  = Wsort(max(1,round(0.05*nv)));
    Whi95 = Wsort(min(nv,round(0.95*nv)));
    Wmean = mean(Wsort);

    fprintf('\n===== Posterior (6 calibration + silicon budget, n = %d samples) =====\n', Npost)
    fprintf('acceptance %.1f%% | W_extra mean %.1f | median %.1f | 5%%-95%% = [%.1f, %.1f] °C·Myr\n', ...
        100*accept/N_iter, Wmean, Wmed, Wlo5, Whi95)
    fprintf('Sign positive over the whole posterior: %d\n', all(Post_W>0))
else
    Wmed=nan; Wlo5=nan; Whi95=nan; Wmean=nan;
    fprintf('\n===== Posterior (6 calibration + silicon budget, n = %d samples; W_extra not computed) =====\n', Npost)
    fprintf('acceptance %.1f%%\n', 100*accept/N_iter)
end
fprintf('Posterior mean theta: k_rw = %.2f Tmol | Q_rw = %.4f | Alk:Si = %.2f\n', ...
    mean(Post(:,1))/1e12, mean(Post(:,2)), mean(Post(:,3)))
fprintf('Reference: prior propagation 19.0 [2.1, 38.6]; pCO2+DSi MCMC 18.6 [2.0, 38.0]\n')

% ---- plotting ----
set(0,'DefaultFigureVisible','on')
if DO_WEXTRA
    figure('Name','6-term calibration + silicon budget posterior','Color','w','Position',[60 60 1100 400])
    subplot(1,3,1); histogram(Post_W,25,'FaceColor',[0.35 0.55 0.85]); hold on
    xline(Wmed,'--r',sprintf('med=%.1f',Wmed)); xline(Wlo5,'--k'); xline(Whi95,'--k')
    xlabel('W_{extra} (°C·Myr)'); ylabel('Sample count'); grid on
    title(sprintf('Posterior W_{extra} [%.1f, %.1f]', Wlo5, Whi95))
    subplot(1,3,2); scatter(Post(:,1)/1e12, Post_W,8,'filled'); xlabel('k_{rw} (Tmol)'); ylabel('W_{extra}'); grid on; title('k_{rw} vs W')
    subplot(1,3,3); scatter(Post(:,2), Post_W,8,'filled'); xlabel('Q_{rw}'); ylabel('W_{extra}'); grid on; title('Q_{rw} vs W')
else
    figure('Name','6-term calibration posterior theta','Color','w','Position',[60 60 700 400])
    subplot(1,2,1); histogram(Post(:,1)/1e12,25); xlabel('k_{rw} (Tmol)'); ylabel('Sample count');
        grid on; title('Posterior k_{rw}')
    subplot(1,2,2); histogram(Post(:,3),25); xlabel('Alk:Si'); ylabel('Sample count'); grid on;
        title('Posterior Alk:Si')
end

save(fullfile(outdir,'calibrate_MCMC.mat'), 'Post','Post_W','Wmed','Wlo5','Whi95','Wmean', ...
    'accept','N_iter','burn','lb','ub','maxSink','t_start','t_end', ...
    'pCO2_obs','DSi_obs','TAlk_obs','T_obs','pH_obs','Om_obs', ...
    'Si_silw_obs','Si_sfw_obs','Si_bio_obs','Si_rw_obs','DO_WEXTRA')
fprintf('\nResults saved to calibrate_MCMC.mat\n')