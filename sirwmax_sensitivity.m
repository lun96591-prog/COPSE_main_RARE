% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Sensitivity of the residual-temperature integral W_extra to Si_rw_max
%%%% Purpose: test how the Si_rw_max cap (a first-order structural sensitivity) affects the clean-window integral.
%%%%       This is the only structural sensitivity that is both first-order and runs cleanly
%%%%       (abiotic and biogenic silica each have caveats).
%%%%
%%%% Logic: Si_rw_max determines whether Early Triassic Si_rw is clamped → controls the DSi substrate
%%%% → controls the temperature-dependent amplification.
%%%% Fixed Q_rw = 0.0608, k_rw = 1.0e12, Alk:Si = 2.66; scan Si_rw_max ∈ {15, 20, 30, 40} Tmol.
%%%% W_extra = ∫(T_Qrw − T_Q0) dt over the clean window (-248.7 → -237 Ma).

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

% ---- parameters ----
Q_ref   = 0.0608;       % physical Q_rw
k_rw    = 1.0e12;       % chassis A
alkrw   = 2.66;
Smax_list = [15, 20, 30, 40] * 1e12;   % Si_rw_max scan settings (mol/yr)
t_start = -248.7;       % integration start (peak)
t_end   = -237;         % clean-window end (Early + Middle Triassic)

nS = numel(Smax_list);
W        = nan(1,nS);
Si_rw_ET = nan(1,nS);   % Early Triassic mean Si_rw (to inspect the cap)
cap_hit  = false(1,nS); % whether the peak reaches the cap

for s = 1:nS
    Smax = Smax_list(s);
    fprintf('Si_rw_max = %.0f Tmol ... ', Smax/1e12)
    try
        % (1) sample
        SENS = struct();
        SENS.k_rw_si        = k_rw;
        SENS.Q_rw           = Q_ref;
        SENS.Alk_Si_RW_ratio = alkrw;
        SENS.Si_rw_max      = Smax;
        run = COPSE_frontend(0);
        close all
        ts  = run.state.time_myr(:);
        Tsm = run.state.tempC(:);
        Srw = run.state.Si_rw(:);

        % (2) paired baseline (Q_rw = 0)
        SENS.Q_rw = 0.0;
        runp = COPSE_frontend(0);
        close all
        tp = runp.state.time_myr(:);
        Tp = runp.state.tempC(:);
        Tp_int = interp1(tp, Tp, ts, 'linear', 'extrap');

        % clean-window residual-temperature integral
        sel = (ts >= t_start) & (ts <= t_end);
        W(s) = trapz(ts(sel), Tsm(sel) - Tp_int(sel));

        % Early Triassic mean Si_rw + cap test
        et = (ts >= -253) & (ts <= -248);
        Si_rw_ET(s) = mean(Srw(et)) / 1e12;
        cap_hit(s) = max(Srw(et)) >= Smax * 0.99;

        fprintf('W=%.1f, ET Si_rw=%.2f (cap=%d)\n', W(s), Si_rw_ET(s), cap_hit(s))
    catch ME
        fprintf('FAIL: %s\n', ME.message)
    end
end

% ---- print ----
fprintf('\n===== Sensitivity of W_extra to Si_rw_max (clean window) =====\n')
fprintf('%-12s %-14s %-14s %-8s\n', 'Si_rw_max', 'W_extra', 'ET Si_rw', 'cap')
for s = 1:nS
    fprintf('%-12.0f %-14.1f %-14.2f %-8d\n', Smax_list(s)/1e12, W(s), Si_rw_ET(s), cap_hit(s))
end
if all(~isnan(W))
    fprintf('\nW_extra range: [%.1f, %.1f] °C·Myr (Si_rw_max 15→40)\n', min(W), max(W))
end

% ---- plotting ----
set(0,'DefaultFigureVisible','on')
figure('Name','Si_rw_max sensitivity','Color','w')
plot(Smax_list/1e12, W, 'o-', 'LineWidth',1.5, 'MarkerSize',8, 'Color',[0 0.4 0.8])
xlabel('Si_{rw,max} (Tmol/yr)'); ylabel('W_{extra} (°C·Myr)')
grid on
title('Sensitivity of the clean-window residual integral to Si_{rw,max}')

save(fullfile(outdir,'sirwmax_sensitivity.mat'), 'Smax_list', 'W', 'Si_rw_ET', 'cap_hit', 't_start', 't_end')
fprintf('\nResults saved to sirwmax_sensitivity.mat\n')