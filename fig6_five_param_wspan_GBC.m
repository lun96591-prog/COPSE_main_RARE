%% ===================================================================
%%  Fig. 6  One-at-a-time sensitivity of W_extra to five parameters
%%          (GBC style)
%% -------------------------------------------------------------------
%%  Changes vs the original five_param_wspan.m:
%%    * in-figure title removed (AGU/GBC: no labels inside the image)
%%    * Arial lettering; Okabe-Ito colour-blind-safe palette
%%    * range annotations kept, but the legend is re-anchored so it no
%%      longer overlaps the range text
%%    * 17.8 cm width, 600 dpi, RGB, white background
%%  The scan definitions and the W_extra integral are unchanged.
%% ===================================================================

clear; close all; clc
clear global SENS Gtune
global SENS Gtune
Gtune = [];
set(0,'DefaultFigureVisible','off')

chassis = 'D:\Program Files\MATLAB\COPSE-main-phase3a-RW temperature-Isson dependence';
outdir  = 'C:\Users\陆宁\Desktop\终极稿件\figures_GBC';
if ~exist(outdir,'dir'), mkdir(outdir); end
cd(chassis)
cachefile = fullfile(outdir,'fig6_five_param_wspan_GBC.mat');

w0 = -248.7;  w1 = -237;      % W_extra window

base = struct();
base.k_rw_si           = 1.0e12;
base.Alk_Si_RW_ratio   = 2.66;
base.Si_rw_max         = 15e12;
base.Q_rw              = 0.0608;
base.Si_abiotic_thresh = 0.8e-3;
base.f_si_evo_PTB      = 0.04;

W_median   = 19.0;            % ensemble median (main-text result)
W_baseline = 10.87;           % calibrated baseline (recomputed below)

%% ---------------------- scan definitions --------------------------
names = { ...
    'Q_{rw} [0.02–0.12] K^{-1}', ...
    'Alk:Si [1.33–4.0]', ...
    'T_{ref} [287–291] K', ...
    'Si_{rw}^{max} [15–40] Tmol/yr', ...
    '\partial [0.6–0.9] mM'};

fields = {'Q_rw','Alk_Si_RW_ratio','T_ref','Si_rw_max','Si_abiotic_thresh'};

values = { ...
    [0.02 0.05 0.0608 0.09 0.12], ...
    [1.33 2.0 2.66 3.0 4.0], ...
    [287 288 289 290 291], ...
    [15 20 30 40]*1e12, ...
    [0.6 0.7 0.8 0.9]*1e-3};

baseline_value = [base.Q_rw, base.Alk_Si_RW_ratio, 289, ...
                  base.Si_rw_max, base.Si_abiotic_thresh];

nS = numel(names);

%% ---------------------- one-dimensional scans ---------------------
W_all  = cell(1,nS);
W_min  = nan(1,nS);
W_max  = nan(1,nS);
W_base = nan(1,nS);

if isfile(cachefile)
    % restyling only: reuse the scan results, no model runs
    fprintf('Loading cached scan results: %s\n', cachefile);
    load(cachefile);
else
for s = 1:nS
    fn   = fields{s};
    vals = values{s};

    shared_base = any(strcmp(fn, {'Q_rw','T_ref'}));

    if shared_base
        SENS = struct();
        SENS.k_rw_si           = base.k_rw_si;
        SENS.Alk_Si_RW_ratio   = base.Alk_Si_RW_ratio;
        SENS.Si_rw_max         = base.Si_rw_max;
        SENS.Si_abiotic_thresh = base.Si_abiotic_thresh;
        SENS.f_si_evo_PTB      = base.f_si_evo_PTB;
        SENS.Q_rw              = 0.0;
        rB = COPSE_frontend(0); close all
    end

    Wv = nan(size(vals));
    for k = 1:numel(vals)
        fprintf('[%d/%d] %s = %.4g ... ', s, nS, fn, vals(k))

        SENS = struct();
        SENS.k_rw_si           = base.k_rw_si;
        SENS.Alk_Si_RW_ratio   = base.Alk_Si_RW_ratio;
        SENS.Si_rw_max         = base.Si_rw_max;
        SENS.Si_abiotic_thresh = base.Si_abiotic_thresh;
        SENS.f_si_evo_PTB      = base.f_si_evo_PTB;
        SENS.Q_rw              = base.Q_rw;
        SENS.(fn)              = vals(k);
        rS = COPSE_frontend(0); close all

        if ~shared_base
            SENS.Q_rw = 0.0;
            rB = COPSE_frontend(0); close all
        end

        tg = (w0:0.02:w1)';
        Ts = interp1(rS.state.time_myr(:), rS.state.tempC(:), tg, 'linear');
        Tb = interp1(rB.state.time_myr(:), rB.state.tempC(:), tg, 'linear');
        Wv(k) = trapz(tg, Ts - Tb);
        fprintf('W_extra = %.2f °C·Myr\n', Wv(k))
    end

    W_all{s}  = Wv;
    W_min(s)  = min(Wv);
    W_max(s)  = max(Wv);
    ib        = find(abs(vals - baseline_value(s)) < 1e-9, 1);
    if ~isempty(ib), W_base(s) = Wv(ib); end
end
W_baseline = mean(W_base(~isnan(W_base)));
save(cachefile, 'names','fields','values','W_all','W_min','W_max','W_base', ...
    'W_baseline','W_median','base','w0','w1');
end

fprintf('\n================ one-at-a-time range summary ================\n');
fprintf('%-28s %-18s %-10s %-10s %-10s\n','parameter','sampled range','W_min','W_max','span');
for s = 1:nS
    fprintf('%-28s %-18s %-10.2f %-10.2f %-10.2f\n', ...
        names{s}, sprintf('%.4g .. %.4g', min(values{s}), max(values{s})), ...
        W_min(s), W_max(s), W_max(s)-W_min(s));
end
fprintf('\nCalibrated W_extra = %.2f °C·Myr | ensemble median = %.1f\n', ...
    W_baseline, W_median);

%% ---------------------- style -------------------------------------
cols = [0.000 0.447 0.741;    % Okabe-Ito blue
        0.835 0.369 0.000;    % vermillion
        0.008 0.620 0.451;    % bluish green
        0.494 0.184 0.556;    % purple
        0.500 0.500 0.500];   % grey

[~, order] = sort(W_max - W_min, 'descend');

figW = 7.0;  figH = 3.4;
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 figW figH]);
set(fig, 'DefaultAxesFontName','Arial', 'DefaultTextFontName','Arial', ...
         'DefaultLegendFontName','Arial', 'DefaultAxesFontSize',8, ...
         'DefaultTextFontSize',8, 'DefaultLegendFontSize',8, ...
         'DefaultAxesLineWidth',0.6, 'DefaultAxesTickDir','in', ...
         'DefaultAxesLayer','top');

ax = axes('Parent',fig,'Position',[0.225 0.165 0.725 0.800]); hold(ax,'on'); box(ax,'on')

hSpan = gobjects(1,nS);
hPts  = gobjects(1,nS);

for row = 1:nS
    s = order(row);
    y = nS + 1 - row;

    hSpan(row) = plot(ax, [W_min(s) W_max(s)], [y y], '-', ...
        'Color', cols(s,:), 'LineWidth', 5.5);
    hPts(row) = plot(ax, W_all{s}, y*ones(size(W_all{s})), 'o', ...
        'MarkerFaceColor','w', 'MarkerEdgeColor', cols(s,:), ...
        'MarkerSize', 4, 'LineWidth', 0.6);
    if ~isnan(W_base(s))
        plot(ax, W_base(s), y, '|', 'Color','k', 'MarkerSize', 8, 'LineWidth', 0.9);
    end
    dlt = W_max(s) - W_min(s);
    if dlt < 0.1
        lab = sprintf('%.3f–%.3f  (Δ%.3f)', W_min(s), W_max(s), dlt);
    else
        lab = sprintf('%.1f–%.1f  (Δ%.1f)', W_min(s), W_max(s), dlt);
    end
    text(ax, W_max(s)+0.55, y, lab, 'FontSize', 8, 'VerticalAlignment','middle')
end

xl = [min(W_min)-3, max(W_max)+14];
xlim(ax, xl)
xline(ax, W_median,   '-', 'Color',[0.60 0.00 0.25], 'LineWidth', 1.0);
xline(ax, W_baseline, ':', 'Color',[0.35 0.35 0.35], 'LineWidth', 0.7);

ylim(ax, [0.4, nS+0.85])
yticks(ax, 1:nS)
yticklabels(ax, flip(names(order)))
ylabel(ax, 'swept parameter (range)')
xlabel(ax, 'W_{extra} (°C·Myr)')
grid(ax,'on')

text(ax, W_median, nS+0.48, sprintf('ensemble median %.1f', W_median), ...
    'Color',[0.60 0.00 0.25], 'FontSize', 8, 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','left')
text(ax, W_baseline, nS+0.48, sprintf('calibrated set %.1f', W_baseline), ...
    'Color',[0.35 0.35 0.35], 'FontSize', 8, 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','right')

legend(ax, [hSpan(1) hPts(1)], {'single-parameter range','sampled values'}, ...
    'Location','southeast', 'FontSize', 8)

%% ---------------------- save --------------------------------------
set(findall(fig,'-property','FontName'), 'FontName', 'Arial');
set(fig,'PaperPositionMode','auto','InvertHardcopy','on','Color','w');
print(fig, fullfile(outdir,'fig6_five_param_wspan_GBC.png'), '-dpng', '-r600');

fprintf('\nSaved: %s\n', fullfile(outdir,'fig6_five_param_wspan_GBC.png'));
