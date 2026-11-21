%%%% Phase 3a | K_sp_cal -> W_extra paired scan
%%%% Purpose: test whether the alkalinity->Omega->mccb channel gain (W_extra)
%%%%          depends on the effective calcite solubility product.
%%%% For each K_sp_cal: two paired runs (Q_rw = 0.0608 vs 0.0608->0), same everything else.
%%%% W_extra = integral over [-248.7, -237] Ma of (T_on - T_off)  [degC*Myr]
%%%% Run from: the project workspace (script cd's into the live model folder)

clear; close all; clc

modelDir = 'D:\Program Files\MATLAB\COPSE-main-phase3a-RW temperature-Isson dependence';
if ~exist(modelDir,'dir'), error('model dir not found: %s', modelDir); end
oldDir = pwd;
cd(modelDir);

global SENS Gtune
Gtune = [];
set(0,'DefaultFigureVisible','off')

Q_on      = 0.0608;
Ksp_list  = [6e-7, 7e-7, 8.4e-7, 1.0e-6, 1.15e-6, 1.3e-6];
win       = [-248.7, -237];      % W_extra window
modern    = [-1, 0];             % modern diagnostic window
ext       = [-253, -248];        % extinction interval

nK = numel(Ksp_list);
% columns: Ksp | pCO2_mod | T_mod | Alk_mod | Om_mod | pCO2_ET | W_extra | meanW | dT_peak | fT_max
R = nan(nK, 10);
ok = false(nK,1);

fprintf('\n=== Phase 3a K_sp -> W_extra paired scan ===\n');
fprintf('Q_on=%.4f | window=[%.1f %.1f] Ma | %d K_sp values x 2 runs\n\n', Q_on, win(1), win(2), nK);

for i = 1:nK
    Ksp = Ksp_list(i);
    fprintf('[%d/%d] K_sp_cal = %.2e ... ', i, nK, Ksp);

    try
        SENS = struct('Q_rw', Q_on, 'K_sp_cal', Ksp);
        run_on = COPSE_frontend(0); close all
        if abs(run_on.pars.K_sp_cal - Ksp)/Ksp > 1e-6
            warning('K_sp_cal override not applied (actual %.3e)', run_on.pars.K_sp_cal);
        end

        SENS = struct('Q_rw', 0, 'K_sp_cal', Ksp);
        run_off = COPSE_frontend(0); close all

        s1 = run_on.state;  s2 = run_off.state;

        m  = (s1.time_myr >= modern(1)) & (s1.time_myr <= modern(2));
        pCO2m = mean(s1.RCO2(m)) * 280;
        Tm    = mean(s1.tempC(m));
        ALKm  = mean(s1.Alk_conc_mM(m));
        Omm   = mean(s1.Omega_cal(m));

        e  = (s1.time_myr >= ext(1)) & (s1.time_myr <= ext(2));
        pCO2e = mean(s1.RCO2(e)) * 280;

        w  = (s1.time_myr >= win(1)) & (s1.time_myr <= win(2));
        t  = s1.time_myr(w);
        T1 = s1.tempC(w);
        T2 = interp1(s2.time_myr, s2.tempC, t, 'linear', 'extrap');
        dT = T1(:) - T2(:);
        W  = trapz(t(:), dT);
        meanW = W / (win(2) - win(1));

        fTm = NaN;
        if isfield(s1,'f_T_rw'), fTm = max(s1.f_T_rw(w)); end

        R(i,:) = [Ksp, pCO2m, Tm, ALKm, Omm, pCO2e, W, meanW, max(dT), fTm];
        ok(i)  = true;
        fprintf('done  W_extra=%6.2f  meanW=%5.2f  pCO2_mod=%4.0f\n', W, meanW, pCO2m);

    catch ME
        fprintf('FAILED: %s\n', ME.message);
    end
end

fprintf('\n%s\n', repmat('-',1,92));
fprintf('%-9s %8s %6s %7s %5s %8s %9s %8s %8s\n', ...
    'K_sp_cal','pCO2mod','Tmod','Alkmod','Ommod','pCO2ET','W_extra','meanW','dTpeak');
fprintf('%s\n', repmat('-',1,92));
for i = 1:nK
    if ok(i)
        fprintf('%-9.2e %8.0f %6.2f %7.3f %5.2f %8.0f %9.2f %8.2f %8.2f\n', R(i,1:9));
    else
        fprintf('%-9.2e   (failed)\n', Ksp_list(i));
    end
end
fprintf('%s\n', repmat('-',1,92));

% ---- save ----
outDir = 'C:\Users\陆宁\AppData\Roaming\TRAE SOLO CN\ModularData\ai-agent\work-mode-projects\6aa37cbc7ee7ee7730e0d042';
T = array2table(R, 'VariableNames', {'K_sp_cal','pCO2_modern','T_modern','Alk_modern','Omega_modern', ...
     'pCO2_extinction','W_extra','mean_warming','dT_peak','f_T_max'});
writetable(T, fullfile(outDir,'Phase3a_Ksp_Wextra_scan_results.csv'));
save(fullfile(outDir,'Phase3a_Ksp_Wextra_scan_results.mat'), 'R','Ksp_list','win','Q_on');
fprintf('saved: Phase3a_Ksp_Wextra_scan_results.csv / .mat\n');

cd(oldDir);
