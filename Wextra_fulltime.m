% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Full-period W_extra curves for Q_rw on/off
%%%%   (a) ΔT(t) = T(on) − T(off)   —— instantaneous integrand of W_extra
%%%%   (b) cumulative integral W(t) = ∫ΔT dt    —— accumulation of W_extra through time
%%%% The x-axis spans the full model run (−900 to 0 Ma), with the ET interval and the W_extra window marked.
%%%% Flat script, no functions; running the whole file with F5 is sufficient.

clear
close all
clc
clear global SENS Gtune

global SENS
global Gtune
Gtune = [];

set(0,'DefaultFigureVisible','off')

chassis = '';   % <-- set this to the path of the COPSE-main-RARE root folder
outdir  = fullfile(chassis,'outputs');
if ~exist(outdir,'dir'), mkdir(outdir); end
if isempty(chassis)
    error('Set chassis to the path of the COPSE-main-RARE root folder before running.');
end
cd(chassis)

% ---- fixed chassis ----
krw0  = 1.0e12;
alk0  = 2.66;
cap0  = 15e12;
th0   = 0.8e-3;
ptb0  = 0.04;
Q_on  = 0.0608;
Q_off = 0;

% ---- paired run: Q_rw on ----
SENS = struct('k_rw_si',krw0,'Alk_Si_RW_ratio',alk0,'Si_rw_max',cap0, ...
    'Si_abiotic_thresh',th0,'f_si_evo_PTB',ptb0,'Q_rw',Q_on);
r  = COPSE_frontend(0); close all
t_on  = r.state.time_myr(:);
T_on  = r.state.tempC(:);
Srw_on = r.state.Si_rw(:)/1e12;
DSi_on = r.state.DSi_conc_uM(:);

% ---- paired run: Q_rw off ----
SENS.Q_rw = Q_off;
r  = COPSE_frontend(0); close all
t_off = r.state.time_myr(:);
T_off = r.state.tempC(:);

% ---- common grid (use the on-run grid and interpolate the off run onto it) ----
T_off_i = interp1(t_off, T_off, t_on, 'linear');
dT   = T_on - T_off_i;                 % instantaneous extra warming (integrand)
Wcum = cumtrapz(t_on, dT);             % cumulative residual-temperature integral

% ---- windows ----
et0 = -253;  et1 = -248;               % ET interval
w0  = -248.7;  w1 = -237;              % W_extra window

iET = (t_on>=et0)&(t_on<=et1);
iW  = (t_on>=w0)&(t_on<=w1);

% ---- key quantities ----
k_dTmax = find(dT==max(dT),1);
fprintf('[1] full-run max ΔT = %.2f °C @ %.1f Ma\n', dT(k_dTmax), t_on(k_dTmax));
fprintf('[2] mean ΔT in the ET interval = %.2f °C; mean ΔT in the W_extra window = %.2f °C\n', ...
    mean(dT(iET)), mean(dT(iW)));
fprintf('[3] full-run cumulative W(t=0) = %.1f °C·Myr\n', Wcum(end));
fprintf('[4] net integral over the W_extra window = %.2f °C·Myr\n', Wcum(find(iW,1,'last')) - Wcum(find(iW,1,'first')));
fprintf('[5] fraction of the full run with ΔT > 0 = %.3f\n', mean(dT>0));

% ================= plotting =================
set(0,'DefaultFigureVisible','on')
fig = figure('Name','Full-period W_extra curves, Q_rw on/off','Color','w','Position',[70 60 980 640]);

% ---- (a) ΔT(t) ----
subplot(2,1,1)
plot(t_on, dT, '-', 'Color', [0.20 0.35 0.65], 'LineWidth', 1.0); hold on
yline(0, '-k');
% window shading
yl = ylim;
patch([et0 et1 et1 et0], [yl(1) yl(1) yl(2) yl(2)], [0.98 0.80 0.50], 'FaceAlpha', 0.25, 'EdgeColor','none');
patch([w0 w1 w1 w0], [yl(1) yl(1) yl(2) yl(2)], [0.80 0.80 0.80], 'FaceAlpha', 0.30, 'EdgeColor','none');
plot(t_on, dT, '-', 'Color', [0.20 0.35 0.65], 'LineWidth', 1.1);
plot(t_on(k_dTmax), dT(k_dTmax), 'o', 'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
text(t_on(k_dTmax), dT(k_dTmax), sprintf('  max %.2f °C @ %.1f Ma', dT(k_dTmax), t_on(k_dTmax)), 'FontSize', 8);
ylabel('\DeltaT = T(Q_{rw}>0) - T(Q_{rw}=0) (°C)')
xlabel('time (Ma)')
title('(a) instantaneous extra warming (integrand of W_{extra})')
legend({'ΔT','0','ET interval','W_{extra} window'}, 'Location','northeast', 'FontSize', 8)
grid on

% ---- (b) cumulative W(t) ----
subplot(2,1,2)
plot(t_on, Wcum, '-', 'Color', [0.60 0.15 0.40], 'LineWidth', 1.2); hold on
yline(0, '-k');
yl = ylim;
patch([et0 et1 et1 et0], [yl(1) yl(1) yl(2) yl(2)], [0.98 0.80 0.50], 'FaceAlpha', 0.25, 'EdgeColor','none');
patch([w0 w1 w1 w0], [yl(1) yl(1) yl(2) yl(2)], [0.80 0.80 0.80], 'FaceAlpha', 0.30, 'EdgeColor','none');
plot(t_on, Wcum, '-', 'Color', [0.60 0.15 0.40], 'LineWidth', 1.3);
ylabel('cumulative \int \DeltaT dt (°C·Myr)')
xlabel('time (Ma)')
title('(b) cumulative residual-temperature integral W(t)')
text(t_on(1)+10, Wcum(end)*0.75, sprintf('W(t=0) = %.1f °C·Myr', Wcum(end)), 'FontSize', 9)
grid on

% ---- save ----
exportgraphics(fig, fullfile(outdir,'Wextra_fulltime.png'), 'Resolution',300);
save(fullfile(outdir,'Wextra_fulltime.mat'), ...
    't_on','T_on','T_off_i','dT','Wcum','Srw_on','DSi_on','et0','et1','w0','w1');
fprintf('\nFigure saved to Wextra_fulltime.png\n');