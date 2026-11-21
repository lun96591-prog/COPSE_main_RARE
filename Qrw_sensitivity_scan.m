% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Q_rw sensitivity scan (Early Triassic response + Ea/Q10 mapping)
%%%% Purpose: scan the reverse-weathering temperature-dependence coefficient Q_rw
%%%%       (corresponding to 14–83 kJ/mol) and quantify its amplification of Early Triassic
%%%%       Si_rw / pCO2 / temperature relative to the Q_rw = 0 baseline.
%%%% Depends on COPSE_frontend(0) from this release.
%%%%
%%%% Conversion relations (GEOCARB convention, reference temperature 288 K):
%%%%   Ea(kJ/mol) = Q_rw * R * 288^2 / 1000    (R=8.314)
%%%%   Q10         = exp(Q_rw * 10)
%%%%   Q_rw = 0.02/0.05/0.0608/0.09/0.12
%%%%        → Ea ≈ 14/34/42/62/83 kJ/mol; Q10 ≈ 1.22/1.65/1.84/2.46/3.32
%%%% Note: f_T_rw is read directly from s.f_T_rw (model T_ref = 289 K); no reference temperature is hard-coded.
%%%%       each setting is a full -1000 → 0 Ma deep-time integration, so six settings take a while.

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

Q_list = [0.0, 0.02, 0.05, 0.0608, 0.09, 0.12];   % Q_rw scan points (K^-1)
win_lo = -253;                       % Early Triassic window start (Ma)
win_hi = -248;                       % Early Triassic window end (Ma)

R_gas  = 8.314;
T_Ea   = 288;                        % Q_rw → Ea reference temperature (GEOCARB convention)

nQ = length(Q_list);
Si_rw_ET = zeros(1,nQ);
pCO2_ET  = zeros(1,nQ);
T_ET     = zeros(1,nQ);
fT_ET     = zeros(1,nQ);
pCO2_mod = zeros(1,nQ);
T_mod    = zeros(1,nQ);
fT_mod   = zeros(1,nQ);

for k = 1:nQ
    Q = Q_list(k);
    fprintf('Running Q_rw = %.4f ...\n', Q)

    SENS = struct();
    SENS.Q_rw = Q;

    run = COPSE_frontend(0);
    close all

    s = run.state;
    idx  = (s.time_myr >= win_lo) & (s.time_myr <= win_hi);
    idxm = (s.time_myr >= -1) & (s.time_myr <= 0);

    Si_rw_ET(k) = mean(s.Si_rw(idx)) / 1e12;    % Tmol Si/yr
    pCO2_ET(k)  = mean(s.RCO2(idx)) * 280;      % ppm
    T_ET(k)     = mean(s.tempC(idx));           % °C
    fT_ET(k)    = mean(s.f_T_rw(idx));          % read directly from the model, not hard-coded
    pCO2_mod(k) = mean(s.RCO2(idxm)) * 280;     % modern self-check (ppm)
    T_mod(k)    = mean(s.tempC(idxm));          % modern self-check (°C)
    fT_mod(k)   = mean(s.f_T_rw(idxm));         % modern f_T_rw self-check (should be ≈1)
end

% increments relative to Q_rw = 0
dT_ET    = T_ET    - T_ET(1);
dpCO2_ET = pCO2_ET - pCO2_ET(1);
dSi_ET   = Si_rw_ET - Si_rw_ET(1);

% Ea / Q10 mapping
Ea_kJ = Q_list * R_gas * T_Ea^2 / 1000;
Q10   = exp(Q_list * 10);

fprintf('\n========== Q_rw sensitivity (with Ea/Q10) ==========\n')
fprintf('Early Triassic window: %d to %d Ma\n', win_lo, win_hi)
fprintf('%-7s %-9s %-7s %-11s %-9s %-6s %-7s %-8s %-8s\n', ...
    'Q_rw','Ea(kJ)','Q10','Si_rw(Tmol)','pCO2(ppm)','T(C)','f_T','dSi_rw','dT(C)')
for k = 1:nQ
    fprintf('%-7.4f %-9.0f %-7.2f %-11.2f %-9.1f %-6.2f %-7.2f %-8.2f %-8.2f\n', ...
        Q_list(k), Ea_kJ(k), Q10(k), Si_rw_ET(k), pCO2_ET(k), T_ET(k), fT_ET(k), dSi_ET(k), dT_ET(k))
end

fprintf('\nModern self-check pCO2 (ppm): '); fprintf('%.1f ', pCO2_mod); fprintf('\n')
fprintf('Modern self-check T    (°C) : '); fprintf('%.2f ', T_mod);    fprintf('\n')
fprintf('Modern self-check f_T_rw     : '); fprintf('%.3f ', fT_mod);    fprintf('\n')

set(0,'DefaultFigureVisible','on')

figure('Name',' Q_rw sensitivity','Color','w')
subplot(1,3,1)
plot(Q_list, Si_rw_ET, 'o-', 'LineWidth', 1.5)
xlabel('Q_{rw} (1/K)'); ylabel('Si_{rw} (Tmol/yr)'); grid on
subplot(1,3,2)
plot(Q_list, dT_ET, 'o-', 'LineWidth', 1.5)
xlabel('Q_{rw} (1/K)'); ylabel('\Delta T (C, vs Q_{rw}=0)'); grid on
subplot(1,3,3)
plot(Q_list, dpCO2_ET, 'o-', 'LineWidth', 1.5)
xlabel('Q_{rw} (1/K)'); ylabel('\Delta pCO_2 (ppm, vs Q_{rw}=0)'); grid on

save(fullfile(outdir,'Qrw_sensitivity.mat'), 'Q_list', 'Ea_kJ', 'Q10', ...
    'Si_rw_ET', 'pCO2_ET', 'T_ET', 'fT_ET', 'dT_ET', 'dpCO2_ET', 'dSi_ET')
fprintf('\nResults saved to Qrw_sensitivity.mat\n')