%% ===================================================================
%%  Fig. 4  Prior propagation vs posterior distributions of W_extra
%%          (GBC style)
%% -------------------------------------------------------------------
%%  Changes vs the original fig2_prior_vs_posterior.m:
%%    * in-figure title removed (AGU/GBC: no labels inside the image)
%%    * Arial lettering; Okabe-Ito colour-blind-safe palette
%%    * reads the Phase3a_* caches that actually exist in the chassis
%%    * 17.8 cm width, 600 dpi, RGB, white background
%%  Plot-only: the model and the MCMC chains are not re-run.
%% ===================================================================

clear; close all; clc

chassis = 'D:\Program Files\MATLAB\COPSE-main-phase3a-RW temperature-Isson dependence';
outdir  = 'C:\Users\陆宁\Desktop\终极稿件\figures_GBC';
if ~exist(outdir,'dir'), mkdir(outdir); end
cd(chassis)

f_prior  = fullfile(chassis, 'Phase3a_warmtail_LHS_calibrated.mat');
f_post2  = fullfile(chassis, 'Phase3a_warmtail_MCMC.mat');
f_post10 = fullfile(chassis, 'Phase3a_calibrate_MCMC.mat');
assert(isfile(f_prior) && isfile(f_post2) && isfile(f_post10), ...
    'Missing result files: %s / %s / %s', f_prior, f_post2, f_post10);

Sp  = load(f_prior);   Wp  = Sp.W(Sp.ok);  Wp = Wp(:);
S2  = load(f_post2);   W2  = S2.Post_W(:);
S10 = load(f_post10);  W10 = S10.Post_W(:);

%% ---------------------- statistics --------------------------------
stat = @(W) [numel(W), median(W), prctile(W,5), prctile(W,95), min(W)];
sp = stat(Wp); s2 = stat(W2); s10 = stat(W10);
fprintf('prior       n=%d median=%.2f 5-95%%=[%.2f, %.2f] min=%.2f\n', sp);
fprintf('posterior2  n=%d median=%.2f 5-95%%=[%.2f, %.2f] min=%.2f\n', s2);
fprintf('posterior10 n=%d median=%.2f 5-95%%=[%.2f, %.2f] min=%.2f\n', s10);

%% ---------------------- style -------------------------------------
c_prior = [0.320 0.480 0.720];
c_p2    = [0.835 0.369 0.000];   % Okabe-Ito vermillion
c_p10   = [0.008 0.620 0.451];   % Okabe-Ito bluish green

figW = 7.0;  figH = 3.3;
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 figW figH]);
set(fig, 'DefaultAxesFontName','Arial', 'DefaultTextFontName','Arial', ...
         'DefaultLegendFontName','Arial', 'DefaultAxesFontSize',8, ...
         'DefaultTextFontSize',8, 'DefaultLegendFontSize',8, ...
         'DefaultAxesLineWidth',0.6, 'DefaultAxesTickDir','in', ...
         'DefaultAxesLayer','top');

ax = axes('Parent',fig,'Position',[0.085 0.175 0.895 0.780]); hold(ax,'on'); box(ax,'on')

xmax  = ceil(max([Wp; W2; W10])/2)*2 + 2;
edges = 0:2:xmax;
dx    = 2;

hP = histogram(ax, Wp, edges, 'Normalization','pdf', 'FaceColor', c_prior, ...
    'EdgeColor', [0.20 0.30 0.50], 'LineWidth', 0.4, 'FaceAlpha', 0.55);

[c2,  ~] = histcounts(W2,  edges);
[c10, ~] = histcounts(W10, edges);
pdf2  = c2  / (numel(W2)  * dx);
pdf10 = c10 / (numel(W10) * dx);

h2  = stairs(ax, edges, [pdf2  pdf2(end)],  '-',  'Color', c_p2,  'LineWidth', 1.2);
h10 = stairs(ax, edges, [pdf10 pdf10(end)], '--', 'Color', c_p10, 'LineWidth', 1.2);

xlim(ax, [0 xmax])
ylabel(ax, 'Probability density')
xlabel(ax, 'W_{extra}  (°C·Myr)')
grid(ax,'on')

legend(ax, [hP h2 h10], { ...
    sprintf('prior propagation, n = %d (median %.1f)',   sp(1),  sp(2)), ...
    sprintf('posterior, 2-term, n = %d (median %.1f)',   s2(1),  s2(2)), ...
    sprintf('posterior, 10-term, n = %d (median %.1f)',  s10(1), s10(2))}, ...
    'Location', 'northeast', 'FontSize', 8);

%% ---------------------- save --------------------------------------
set(findall(fig,'-property','FontName'), 'FontName', 'Arial');
set(fig,'PaperPositionMode','auto','InvertHardcopy','on','Color','w');
print(fig, fullfile(outdir,'fig4_prior_vs_posterior_GBC.png'), '-dpng', '-r600');
save(fullfile(outdir,'fig4_prior_vs_posterior_GBC.mat'), 'Wp','W2','W10','sp','s2','s10');

fprintf('\nSaved: %s\n', fullfile(outdir,'fig4_prior_vs_posterior_GBC.png'));
