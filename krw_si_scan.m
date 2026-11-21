% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% k_rw_si sensitivity scan (case B/C comparison)
%%%% Purpose: scan the reverse-weathering rate constant k_rw_si (1.0 / 1.5 / 2.0 Tmol Si/yr equivalent),
%%%%       comparing the modern steady state (DSi / Si_rw / pCO2 / T / Omega / TAlk)
%%%%       with the Early Triassic window (Si_rw / pCO2 / T / dT vs the chassis),
%%%%       to quantify the climate impact of the literature range of the reverse-weathering flux
%%%%       (estuarine lower bound vs Rahman effective value).
%%%%
%%%% Three settings:
%%%%   A: k_rw_si = 1.0e12  ← current chassis (estuarine value; Tréguer & De La Rocha 2013)
%%%%   B: k_rw_si = 1.5e12  ← Tréguer et al. (2021) synthesis central value (recommended box-model mean)
%%%%   C: k_rw_si = 2.0e12  ← upper effective value of Rahman et al. (2017) after excluding bSi→clay
%%%%
%%%% Key literature background:
%%%%   Tréguer & De La Rocha (2013): estuarine mass balance = 1.5 ± 0.5 Tmol/yr (water-column DSi removal only)
%%%%   Rahman et al. (2017) GBC: 4.5–4.9 Tmol/yr (including 2.7 Tmol of bSi→clay conversion,
%%%%       which the Si_bio sink already covers and which a box model must exclude)
%%%%   Tréguer et al. (2021), Biogeosciences: synthesis = ~2 Tmol/yr (true whole-ocean extra DSi removal)
%%%%
%%%% Dependency: COPSE_frontend(0) from this release.
%%%% Note: k_rw_si is already overridden in the frontend through the global SENS;
%%%%       Si_rw_max = 15e12 is the live chassis value; the scan checks whether the Early Triassic hits the cap.

clear
close all
clc
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
Gtune = [];                          % use the default initial state

set(0,'DefaultFigureVisible','off')  % suppress COPSE_plot pop-ups

% ---- scan settings ----
krw_list = [1.0e12, 1.5e12, 2.0e12];   % A (chassis) / B (recommended) / C (upper bound)  — mol Si/yr
win_lo   = -253;                       % Early Triassic window start (Ma)
win_hi   = -248;                       % Early Triassic window end (Ma)

nK = length(krw_list);

% ---- preallocate result arrays ----
% modern steady state
DSi_mod     = zeros(1,nK);            % μM
Si_rw_mod   = zeros(1,nK);            % Tmol Si/yr
pCO2_mod    = zeros(1,nK);            % ppm
T_mod       = zeros(1,nK);            % °C
Omega_mod   = zeros(1,nK);            % calcite saturation
TAlk_mod    = zeros(1,nK);            % mmol/kg
DIC_mod     = zeros(1,nK);            % mmol/kg
fT_mod      = zeros(1,nK);            % f_T_rw (should be ≈1.0)

% Early Triassic window means
DSi_ET      = zeros(1,nK);            % μM
Si_rw_ET    = zeros(1,nK);            % Tmol Si/yr
pCO2_ET     = zeros(1,nK);            % ppm
T_ET        = zeros(1,nK);            % °C
fT_ET       = zeros(1,nK);            % f_T_rw
cap_hit     = zeros(1,nK);            % whether Early Triassic Si_rw hits the Si_rw_max cap (boolean)

% ---- run ----
for k = 1:nK
    krw = krw_list(k);
    fprintf('Running k_rw_si = %.1e (case %c) ...\n', krw, char('A'+k-1))

    SENS = struct();
    SENS.k_rw_si  = krw;
    SENS.Si_rw_max = 15e12;   % uniform nominal value of 15 Tmol (numerically stable cap)

    run = COPSE_frontend(0);
    close all

    s    = run.state;
    pars = run.pars;

    idx  = (s.time_myr >= win_lo) & (s.time_myr <= win_hi);
    idxm = (s.time_myr >= -1) & (s.time_myr <= 0);

    % --- modern ---
    DSi_mod(k)   = mean(s.DSi_conc_uM(idxm));                   % already in μM
    Si_rw_mod(k) = mean(s.Si_rw(idxm)) / 1e12;                  % mol/yr -> Tmol/yr
    pCO2_mod(k)  = mean(s.RCO2(idxm)) * 280;                    % ratio × 280 ppm
    T_mod(k)     = mean(s.tempC(idxm));                         % °C
    Omega_mod(k) = mean(s.Omega_cal(idxm));
    TAlk_mod(k)  = mean(s.Alk_conc_mM(idxm));                   % already in mmol/kg
    DIC_mod(k)   = mean(s.DIC_conc_mM(idxm));                   % already in mmol/kg
    fT_mod(k)    = mean(s.f_T_rw(idxm));

    % --- Early Triassic ---
    DSi_ET(k)    = mean(s.DSi_conc_uM(idx));
    Si_rw_ET(k)  = mean(s.Si_rw(idx)) / 1e12;                   % mol/yr -> Tmol/yr
    pCO2_ET(k)   = mean(s.RCO2(idx)) * 280;
    T_ET(k)      = mean(s.tempC(idx));
    fT_ET(k)     = mean(s.f_T_rw(idx));

    % --- whether the Si_rw_max cap is reached ---
    cap_hit(k)    = max(s.Si_rw(idx)) >= (pars.Si_rw_max * 0.99);   % 99% threshold
end

% ---- increments relative to A (chassis) ----
dT_ET    = T_ET    - T_ET(1);
dpCO2_ET = pCO2_ET - pCO2_ET(1);
dSi_ET   = Si_rw_ET - Si_rw_ET(1);

krw_label = {'A:1.0','B:1.5','C:2.0'};   % Tmol/yr

% ---- print results ----
fprintf('\n========== k_rw_si sensitivity ==========\n')
fprintf('Early Triassic window: %d to %d Ma\n', win_lo, win_hi)
fprintf('Si_rw_max = %.0f Tmol/yr (live chassis)\n\n', pars.Si_rw_max/1e12)

fprintf('--- modern steady state ---\n')
fprintf('%-7s %-8s %-8s %-9s %-6s %-7s %-7s %-7s %-7s\n', ...
    'case','k_rw(T)','DSi(uM)','Si_rw(T)','pCO2','T(C)','Omega','TAlk','f_T')
for k = 1:nK
    fprintf('%-7s %-8.1f %-8.1f %-9.2f %-6.1f %-7.2f %-7.2f %-7.3f %-7.3f\n', ...
        krw_label{k}, krw_list(k)/1e12, DSi_mod(k), Si_rw_mod(k), ...
        pCO2_mod(k), T_mod(k), Omega_mod(k), TAlk_mod(k), fT_mod(k))
end

fprintf('\n--- Early Triassic window ---\n')
fprintf('%-7s %-8s %-8s %-9s %-7s %-6s %-7s %-7s %-7s\n', ...
    'case','k_rw(T)','DSi(uM)','Si_rw(T)','pCO2','T(C)','dT(C)','dpCO2','cap?')
for k = 1:nK
    fprintf('%-7s %-8.1f %-8.1f %-9.2f %-7.1f %-6.2f %-7.2f %-7.0f %-7d\n', ...
        krw_label{k}, krw_list(k)/1e12, DSi_ET(k), Si_rw_ET(k), ...
        pCO2_ET(k), T_ET(k), dT_ET(k), dpCO2_ET(k), cap_hit(k))
end

% ---- plotting ----
set(0,'DefaultFigureVisible','on')

figure('Name',' k_rw_si sensitivity','Color','w','Position',[50 50 1100 500])

subplot(2,4,1)
plot(krw_list/1e12, DSi_mod, 'o-', 'LineWidth',1.5, 'MarkerSize',8)
xlabel('k_{rw} (Tmol Si/yr)'); ylabel('modern DSi (\muM)'); grid on

subplot(2,4,2)
plot(krw_list/1e12, Si_rw_mod, 'o-', 'LineWidth',1.5, 'MarkerSize',8)
xlabel('k_{rw} (Tmol Si/yr)'); ylabel('modern Si_{rw} (Tmol/yr)'); grid on

subplot(2,4,3)
plot(krw_list/1e12, pCO2_mod, 'o-', 'LineWidth',1.5, 'MarkerSize',8)
xlabel('k_{rw} (Tmol Si/yr)'); ylabel('modern pCO_2 (ppm)'); grid on

subplot(2,4,4)
plot(krw_list/1e12, T_mod, 'o-', 'LineWidth',1.5, 'MarkerSize',8)
xlabel('k_{rw} (Tmol Si/yr)'); ylabel('modern T (°C)'); grid on

subplot(2,4,5)
plot(krw_list/1e12, Omega_mod, 'o-', 'LineWidth',1.5, 'MarkerSize',8, 'Color',[0 0.4 0.8])
hold on; yline(2.0, '--r', '\Omega=2'); hold off
xlabel('k_{rw} (Tmol Si/yr)'); ylabel('modern \Omega_{cal}'); grid on

subplot(2,4,6)
plot(krw_list/1e12, TAlk_mod, 'o-', 'LineWidth',1.5, 'MarkerSize',8, 'Color',[0.8 0.4 0])
hold on; yline(2.3, '--r', '2.3 mM'); hold off
xlabel('k_{rw} (Tmol Si/yr)'); ylabel('modern TAlk (mmol/kg)'); grid on

subplot(2,4,7)
plot(krw_list/1e12, Si_rw_ET, 's-', 'LineWidth',1.5, 'MarkerSize',8, 'Color',[0.6 0 0.8])
hold on; yline(pars.Si_rw_max/1e12, '--k', sprintf('Si_{rw,max}=%.0f',pars.Si_rw_max/1e12)); hold off
xlabel('k_{rw} (Tmol Si/yr)'); ylabel('Early Triassic Si_{rw} (Tmol/yr)'); grid on

subplot(2,4,8)
plot(krw_list/1e12, dT_ET, 's-', 'LineWidth',1.5, 'MarkerSize',8, 'Color',[0.8 0 0])
xlabel('k_{rw} (Tmol Si/yr)'); ylabel('Early Triassic \Delta T vs A (°C)'); grid on

sgtitle('Reverse-weathering strength k_{rw} sensitivity (Isson threshold + Q_{rw}=0.0608 + PTB collapse valley 0.04)')

% ---- save ----
save(fullfile(outdir,'krw_si_scan.mat'), 'krw_list', 'krw_label', ...
    'DSi_mod','Si_rw_mod','pCO2_mod','T_mod','Omega_mod','TAlk_mod','DIC_mod','fT_mod', ...
    'DSi_ET','Si_rw_ET','pCO2_ET','T_ET','fT_ET', ...
    'dT_ET','dpCO2_ET','dSi_ET','cap_hit', 'win_lo','win_hi')

fprintf('\nResults saved to krw_si_scan.mat\n')
