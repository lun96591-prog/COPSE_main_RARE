% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Q_rw sign-symmetry control
%%%%
%%%% Purpose: answer whether W_extra > 0 is merely an artefact of assuming a positive activation energy.
%%%%
%%%% Logical premise (essential to understand, otherwise the run is wasted):
%%%%   under the paired design, the two runs at Q_rw = 0 are **identical**, so W_extra ≡ 0;
%%%%   the zero crossing is an identity, not a result, and needs no distributional "verification".
%%%%   what must be tested is whether the **magnitude is antisymmetric about Q_rw**:
%%%%     if |W(+Q)| ≈ |W(-Q)|, the sign comes purely from the input sign and the **magnitude is the model's own dynamical answer**;
%%%%     if clearly asymmetric, the model's response to warming/cooling is itself asymmetric (a stronger structural finding).
%%%%
%%%% Design: deterministic single-parameter scan rather than a 2000-sample ensemble.
%%%%   The sign question has one degree of freedom; 9 Q_rw values + 1 shared baseline = 10 deep-time integrations ≈ 30 s.
%%%%
%%%% Fixed chassis: k_rw = 1.0 Tmol/yr, Alk:Si = 2.66, Si_rw_max = 15 Tmol/yr,
%%%%   clean window -248.7 → -237 Ma, identical to the main result (only Q_rw varies).

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

% ---- chassis and scan settings ----
krw_base   = 1.0e12;          % calibrated chassis k_rw
alk_base   = 2.66;            % default Alk:Si (sepiolite/smectite)
Srw_max    = 15e12;
t_start    = -248.7;
t_end      = -237;
Tref       = 289;             % consistent with equations
Rgas       = 8.314;

Qlist = [-0.12, -0.09, -0.0608, -0.03, 0.0, 0.03, 0.0608, 0.09, 0.12];  % symmetric range
nQ    = numel(Qlist);

% ---- (1) run the shared baseline first (Q_rw = 0, chassis k_rw/Alk:Si) ----
fprintf('Baseline run Q_rw = 0 ...\n')
SENS = struct();
SENS.k_rw_si         = krw_base;
SENS.Q_rw            = 0.0;
SENS.Alk_Si_RW_ratio = alk_base;
SENS.Si_rw_max       = Srw_max;
run0 = COPSE_frontend(0);
close all
t0 = run0.state.time_myr(:);
T0 = run0.state.tempC(:);

% ---- (2) run a sample for each Q_rw and difference against the baseline ----
Wq = nan(1, nQ);
for q = 1:nQ
    Q = Qlist(q);
    fprintf('[%d/%d] Q_rw = %+.4f ... ', q, nQ, Q)

    if Q == 0
        Wq(q) = 0.0;                      % identity: identical to the baseline
        fprintf('W = 0 (identical runs by construction)\n')
        continue
    end

    SENS = struct();
    SENS.k_rw_si         = krw_base;
    SENS.Q_rw            = Q;
    SENS.Alk_Si_RW_ratio = alk_base;
    SENS.Si_rw_max       = Srw_max;
    run = COPSE_frontend(0);
    close all

    ts  = run.state.time_myr(:);
    Tsm = run.state.tempC(:);
    T0i = interp1(t0, T0, ts, 'linear', 'extrap');   % interpolate the baseline onto the same grid

    sel = (ts >= t_start) & (ts <= t_end);
    Wq(q) = trapz(ts(sel), Tsm(sel) - T0i(sel));
    fprintf('W = %+.1f °C·Myr\n', Wq(q))
end

% ---- (3) report: Q_rw / Ea / Q10 / W_extra ----
fprintf('\n===== Q_rw sign symmetry (k_rw=%.1f Tmol, Alk:Si=%.2f, cap=%.0f, window -248.7→-237) =====\n', ...
    krw_base/1e12, alk_base, Srw_max/1e12)
fprintf('%-10s %-10s %-8s %-12s\n', 'Q_rw(1/K)', 'Ea(kJ/mol)', 'Q10', 'W_extra')
fprintf('%s\n', repmat('-', 1, 46))
for q = 1:nQ
    Q = Qlist(q);
    fprintf('%+10.4f %+10.1f %8.2f %+12.1f\n', ...
        Q, Q*Rgas*Tref^2/1000, exp(10*Q), Wq(q))
end

% ---- (4) antisymmetry test: |W(+Q)| vs |W(-Q)| ----
fprintf('\nAntisymmetry test (pairs with equal |Q|):\n')
fprintf('%-10s %-12s %-12s %-10s\n', '|Q_rw|', 'W(+Q)', 'W(-Q)', 'W(+)+W(-)')
for Q = unique(abs(Qlist(Qlist ~= 0)))
    ip = find(abs(Qlist - Q)     < 1e-12);
    im = find(abs(Qlist + Q)     < 1e-12);
    if ~isempty(ip) && ~isempty(im)
        fprintf('%-10.4f %+12.1f %+12.1f %+10.2f\n', Q, Wq(ip), Wq(im), Wq(ip)+Wq(im))
    end
end
fprintf('\nInterpretation: if W(+) + W(-) is much smaller than |W|, the response is approximately antisymmetric ——\n')
fprintf('      the positive sign comes from the assumed input sign; the magnitude comes from the model dynamics.\n')
if all(Wq(Qlist>0) > 0) && all(Wq(Qlist<0) < 0)
    fprintf('Direction: Q_rw > 0 → W_extra all positive; Q_rw < 0 → all negative (sign flips with the input, as expected).\n')
else
    fprintf('WARNING: samples violate the expected sign; check case by case.\n')
end

% ---- (5) plotting ----
set(0,'DefaultFigureVisible','on')
figure('Name','Q_{rw} sign symmetry','Color','w','Position',[60 60 900 380])
plot(Qlist, Wq, 'o-', 'LineWidth', 1.8, 'MarkerSize', 8, 'Color', [0.2 0.4 0.8])
hold on
yline(0, '-k', 'LineWidth', 1.0)
xline(0, '--k')
xlabel('Q_{rw} (K^{-1})'); ylabel('W_{extra} (°C·Myr)'); grid on
title(sprintf('W_{extra} vs Q_{rw} sign (k_{rw}=%.1f Tmol, Alk:Si=%.2f)', krw_base/1e12, alk_base))

save(fullfile(outdir,'Qrw_sign_symmetry.mat'), 'Qlist', 'Wq', 'krw_base', 'alk_base', ...
    'Srw_max', 't_start', 't_end', 'Tref')
fprintf('\nResults saved to Qrw_sign_symmetry.mat\n')
