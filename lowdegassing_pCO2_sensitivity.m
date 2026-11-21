% NOTE: set `chassis` at the top of this script to the path of your local
% COPSE-main-RARE root folder (the folder containing COPSE_frontend.m and
% COPSE_equations.m) before running.
% Part of the COPSE-main-RARE release (modifications (c) 2026 Ning Lu, Xiamen University).
% Distributed under the Apache License, Version 2.0; see ../LICENSE.
%%%% Low-degassing scenario (flat script, no functions; run entirely with F5 or run())
%%%% Entry point COPSE_frontend(0); scaling via SENS.DEGASS_scale; modern state at t = 0 instantaneous

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
kbio0 = 3.5e12;
Q_on  = 0.0608;
Q_off = 0;

w0 = -248.7;  w1 = -237;  Lwin = w1 - w0;
pC_target = 280;
tg = (w0:0.02:w1)';        % common integration grid

% ---- generic: assemble SENS (rewritten every time to avoid stale fields) ----
SENS = struct();
SENS.k_rw_si         = krw0;
SENS.Alk_Si_RW_ratio = alk0;
SENS.Si_rw_max       = cap0;
SENS.Si_abiotic_thresh = th0;
SENS.f_si_evo_PTB    = ptb0;
SENS.k_bio_si        = kbio0;

% ================= (1) baseline f = 1, Q_rw on =================
SENS.Q_rw = Q_on;
r  = COPSE_frontend(0); close all
t0_on  = r.state.time_myr(:);
T0_on  = r.state.tempC(:);
pC0    = r.state.RCO2(:)*280;
Ds0    = r.state.DSi_conc_uM(:);
[~,j0] = min(abs(t0_on));
iET0 = (t0_on>=-253)&(t0_on<=-248);
pC0_t0 = pC0(j0);  T0_t0 = T0_on(j0);  DSi0_t0 = Ds0(j0);
ET_T0 = mean(T0_on(iET0));  ET_pC0 = mean(pC0(iET0));

% ---- baseline Q_rw off ----
SENS.Q_rw = Q_off;
r  = COPSE_frontend(0); close all
T0_off = r.state.tempC(:);
T0_off_g = interp1(r.state.time_myr(:), T0_off, tg, 'linear');
W0 = trapz(tg, interp1(t0_on,T0_on,tg,'linear') - T0_off_g);

fprintf('BASELINE|f=1|pC_t0=%.1f|T_t0=%.2f|DSi_t0=%.1f|ET_T=%.2f|ET_pC=%.0f|W=%.2f|meanW=%.2f\n', ...
    pC0_t0, T0_t0, DSi0_t0, ET_T0, ET_pC0, W0, W0/Lwin);

% ================= (2) bisect f*: t = 0 instantaneous pCO2 -> 280 =================
f_lo = 0.30;  f_hi = 1.00;
for it = 1:12
    fm = 0.5*(f_lo + f_hi);
    SENS.Q_rw = Q_on;
    SENS.DEGASS_scale = fm;
    r  = COPSE_frontend(0); close all
    tm = r.state.time_myr(:);  pm = r.state.RCO2(:)*280;
    [~,jm] = min(abs(tm));
    pC_m = pm(jm);
    fprintf('  bisect f=%.4f -> pC(t0)=%.1f ppm\n', fm, pC_m);
    if pC_m > pC_target, f_hi = fm; else, f_lo = fm; end
end
fstar = 0.5*(f_lo + f_hi);

% ================= (3) paired W at f* (k_bio fixed, no regression) =================
% Note: degassing falls only ~5% and DSi drifts ~2 µM, so regressing k_bio would add projection
% noise; k_bio is therefore fixed directly
SENS.Q_rw = Q_on;
SENS.k_bio_si = kbio0;
SENS.DEGASS_scale = fstar;
r  = COPSE_frontend(0); close all
t2_on = r.state.time_myr(:);  T2_on = r.state.tempC(:);
pC2 = r.state.RCO2(:)*280;  Ds2 = r.state.DSi_conc_uM(:);
[~,j2] = min(abs(t2_on));
pC_c_t0 = pC2(j2);  T_c_t0 = T2_on(j2);  DSi_c_t0 = Ds2(j2);
iET2 = (t2_on>=-253)&(t2_on<=-248);
ET_T2 = mean(T2_on(iET2));  ET_pC2 = mean(pC2(iET2));

SENS.Q_rw = Q_off;
r  = COPSE_frontend(0); close all
T2_off = r.state.tempC(:);
T2_off_g = interp1(r.state.time_myr(:), T2_off, tg, 'linear');
W1 = trapz(tg, interp1(t2_on,T2_on,tg,'linear') - T2_off_g);

% ================= summary =================
fprintf('\n================ Low-degassing scenario (t = 0 instantaneous convention) ================\n');
fprintf('%-24s %-8s %-8s %-8s %-8s %-9s %-9s %-11s\n', ...
    'case','f DEGASS','pCO2_t0','T_t0','DSi_t0','ET_Tmean','ET_pC','W_extra');
fprintf('%-24s %-8.3f %-8.1f %-8.2f %-8.1f %-9.2f %-9.0f %-11.2f\n', ...
    'baseline (pCO2~321)', 1.0, pC0_t0, T0_t0, DSi0_t0, ET_T0, ET_pC0, W0);
fprintf('%-24s %-8.3f %-8.1f %-8.2f %-8.1f %-9.2f %-9.0f %-11.2f\n', ...
    'low-degassing (k_bio fixed)', fstar, pC_c_t0, T_c_t0, DSi_c_t0, ET_T2, ET_pC2, W1);

fprintf('\n[check 1] modern pCO2(t0): %.1f -> %.1f ppm (target %d)\n', pC0_t0, pC_c_t0, pC_target);
fprintf('[check 2] W_extra: %.2f -> %.2f °C·Myr (change %+.1f%%)\n', W0, W1, 100*(W1-W0)/W0);
fprintf('[check 3] mean warming: %.2f -> %.2f °C (window %.1f Myr)\n', W0/Lwin, W1/Lwin, Lwin);
fprintf('[check 4] Early Triassic mean temperature: %.2f -> %.2f °C; pCO2 %.0f -> %.0f ppm\n', ...
    ET_T0, ET_T2, ET_pC0, ET_pC2);
fprintf('[check 5] f* = %.3f (degassing reduced to %.0f%%), k_bio fixed at %.2f e12\n', ...
    fstar, 100*fstar, kbio0/1e12);

save(fullfile(outdir,'lowdegassing_pCO2.mat'), ...
    'fstar','kbio0','W0','W1','pC0_t0','T0_t0','DSi0_t0', ...
    'pC_c_t0','T_c_t0','DSi_c_t0','ET_T0','ET_T2','ET_pC0','ET_pC2', ...
    'w0','w1','pC_target');
fprintf('\nResults saved to lowdegassing_pCO2.mat\n');