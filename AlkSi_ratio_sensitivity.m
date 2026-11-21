% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Alk:Si ratio (Alk_Si_RW_ratio) sensitivity scan
%%%% Purpose: quantify how the stoichiometric coefficient (alkalinity consumed per mol Si by
%%%%       reverse weathering) affects the modern calibration, the Early Triassic response
%%%%
%%%%       and the clean-window residual-temperature integral W_extra. Background: the reverse-
%%%%       weathering alkalinity sink = Alk_Si_RW_ratio × Si_rw, so Alk:Si and k_rw are parallel
%%%%       linear multipliers on that sink and first-order controls on the climate effect; Alk:Si
%%%%       was previously sampled only as one of three LHS parameters. This script scans it explicitly.
%%%%
%%%% Literature range (Zhao et al. 2024 Sci Adv; Rauzi 2024 SM):
%%%%   sepiolite 1.33 → saponite 2.33 → smectite 2.66/2.67(KT2020) → berthierine/chlorite 4.0
%%%%
%%%% Fixed: k_rw = 1.0e12, Q_rw = 0.0608, Si_rw_max = 15e12.
%%%% Each setting runs a sample (Q = 0.0608) plus a paired baseline (Q = 0) for clean-window W_extra.

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

% ---- scan settings ----
ratio_list = [1.33, 2.0, 2.66, 3.0, 4.0];   % Alk:Si (eq Alk per mol Si)
k_rw    = 1.0e12;
Q_ref   = 0.0608;
Srw_max = 15e12;
t_start = -248.7;    % clean-window start (Early Triassic peak)
t_end   = -237;      % clean-window end (avoiding the Late Triassic artifact)

nR = numel(ratio_list);

% preallocation
W_extra  = nan(1,nR);
TAlk_mod = nan(1,nR);  DIC_mod  = nan(1,nR);  Omega_mod = nan(1,nR);
mccb_mod = nan(1,nR);  pCO2_mod = nan(1,nR);  T_mod     = nan(1,nR);
Si_rw_mod= nan(1,nR);  DSi_mod  = nan(1,nR);
TAlk_ET  = nan(1,nR);  pCO2_ET  = nan(1,nR);  T_ET      = nan(1,nR);
Si_rw_ET = nan(1,nR);
Ttraj = cell(1,nR);   % temperature trajectories of the sample runs (Q_rw = 0.0608)
Tbase = cell(1,nR);   % temperature trajectories of the paired baselines (Q_rw = 0, same grid)
tgrid = cell(1,nR);   % time grid

for r = 1:nR
    alkrw = ratio_list(r);
    fprintf('Alk:Si = %.2f ... ', alkrw)
    try
        % (1) sample: Q_rw = Q_ref
        SENS = struct();
        SENS.k_rw_si        = k_rw;
        SENS.Q_rw           = Q_ref;
        SENS.Alk_Si_RW_ratio = alkrw;
        SENS.Si_rw_max      = Srw_max;
        run = COPSE_frontend(0);
        close all
        st  = run.state;
        ts  = st.time_myr(:);
        Tsm = st.tempC(:);

        idxm = (ts >= -1) & (ts <= 0);
        idx  = (ts >= -253) & (ts <= -248);
        TAlk_mod(r) = mean(st.Alk_conc_mM(idxm));
        DIC_mod(r)  = mean(st.DIC_conc_mM(idxm));
        Omega_mod(r)= mean(st.Omega_cal(idxm));
        mccb_mod(r) = mean(st.mccb(idxm))/1e12;
        pCO2_mod(r) = mean(st.RCO2(idxm))*280;
        T_mod(r)    = mean(Tsm(idxm));
        Si_rw_mod(r)= mean(st.Si_rw(idxm))/1e12;
        DSi_mod(r)  = mean(st.DSi_conc_uM(idxm));
        TAlk_ET(r)  = mean(st.Alk_conc_mM(idx));
        pCO2_ET(r)  = mean(st.RCO2(idx))*280;
        T_ET(r)     = mean(Tsm(idx));
        Si_rw_ET(r) = mean(st.Si_rw(idx))/1e12;

        % (2) paired baseline: Q_rw = 0 (same Alk:Si)
        SENS.Q_rw = 0.0;
        runp = COPSE_frontend(0);
        close all
        tp = runp.state.time_myr(:);
        Tp = runp.state.tempC(:);
        Tp_int = interp1(tp, Tp, ts, 'linear', 'extrap');

        % store trajectories for comparison
        tgrid{r} = ts;  Ttraj{r} = Tsm;  Tbase{r} = Tp_int;

        % clean-window residual-temperature integral
        sel = (ts >= t_start) & (ts <= t_end);
        W_extra(r) = trapz(ts(sel), Tsm(sel) - Tp_int(sel));

        fprintf('W=%.1f, mod pCO2=%.0f, ET T=%.2f\n', ...
            W_extra(r), pCO2_mod(r), T_ET(r))
    catch ME
        fprintf('FAIL: %s\n', ME.message)
    end
end

% ---- print ----
fprintf('\n===== Alk:Si ratio sensitivity (k_rw = 1.0, Q_rw = 0.0608, cap = 15) =====\n')
fprintf('%-8s %-11s %-11s %-9s %-8s %-8s %-8s %-8s\n', ...
    'Alk:Si','W_extra','pCO2m','T_mod','TAlk_m','Omega','Si_rw_m','T_ET')
for r = 1:nR
    fprintf('%-8.2f %-11.1f %-11.0f %-9.2f %-8.2f %-8.2f %-8.2f %-8.2f\n', ...
        ratio_list(r), W_extra(r), pCO2_mod(r), T_mod(r), TAlk_mod(r), ...
        Omega_mod(r), Si_rw_mod(r), T_ET(r))
end

fprintf('\nModern-calibration deviation check (relative to 2.66):\n')
i266 = find(ratio_list==2.66);
for r = 1:nR
    fprintf('  %.2f:ΔpCO2 = %+.0f ppm, ΔTAlk = %+.3f mmol/kg, ΔT = %+.2f °C\n', ...
        ratio_list(r), pCO2_mod(r)-pCO2_mod(i266), ...
        TAlk_mod(r)-TAlk_mod(i266), T_mod(r)-T_mod(i266))
end

fprintf('\nW_extra range [%.1f, %.1f]; all settings positive: %s\n', ...
    min(W_extra), max(W_extra), string(all(W_extra>0)))

% ---- plotting ----
set(0,'DefaultFigureVisible','on')
figure('Name','Alk:Si ratio sensitivity','Color','w','Position',[50 50 1100 420])
subplot(1,2,1)
plot(ratio_list, W_extra, 'o-', 'LineWidth',1.5, 'MarkerSize',9, 'Color',[0 0.4 0.8])
xlabel('Alk:Si ratio (eq Alk / mol Si)'); ylabel('W_{extra} (°C·Myr)')
grid on; title('Alk:Si vs clean-window residual-temperature integral')
subplot(1,2,2)
yyaxis left
plot(ratio_list, pCO2_mod, 's-', 'LineWidth',1.5, 'Color',[0.8 0.2 0.2])
ylabel('modern pCO_2 (ppm)')
yyaxis right
plot(ratio_list, T_mod, '^-', 'LineWidth',1.5, 'Color',[0.8 0.5 0])
ylabel('modern T (°C)')
xlabel('Alk:Si ratio'); grid on; title('Sensitivity of the modern calibration to Alk:Si')

% ---- temperature-trajectory comparison (source of the non-monotonic W_extra) ----
pick = [1, 3, 5];            % the 1.33 / 2.66 / 4.0 settings
figure('Name','Alk:Si temperature-trajectory comparison','Color','w','Position',[80 80 1100 500])
for ii = 1:numel(pick)
    r  = pick(ii);
    ts = tgrid{r};  Ts = Ttraj{r};  Tb = Tbase{r};
    sel = (ts >= t_start) & (ts <= t_end);
    subplot(1,3,ii); hold on
    fill([ts(sel); flipud(ts(sel))], [Ts(sel); flipud(Tb(sel))], ...
        [0.60 0.70 0.90], 'EdgeColor','none', 'FaceAlpha',0.35)
    plot(ts, Ts, 'r-', 'LineWidth',1.6)     % Q_rw = 0.0608
    plot(ts, Tb, 'b-', 'LineWidth',1.6)     % Q_rw = 0
    xlim([-250 -236]); xlabel('t (Ma)'); ylabel('T (°C)'); grid on
    title(sprintf('Alk:Si = %.2f   W = %.1f', ratio_list(r), W_extra(r)))
end
legend({'W_{extra} (shaded)','Q_{rw}=0.0608','Q_{rw}=0'}, ...
    'Location','southoutside','Orientation','horizontal')

save(fullfile(outdir,'AlkSi_ratio_sensitivity.mat'), 'ratio_list', 'W_extra', ...
    'TAlk_mod','DIC_mod','Omega_mod','mccb_mod','pCO2_mod','T_mod','Si_rw_mod','DSi_mod', ...
    'TAlk_ET','pCO2_ET','T_ET','Si_rw_ET', 't_start','t_end')
fprintf('\nResults saved to AlkSi_ratio_sensitivity.mat\n')