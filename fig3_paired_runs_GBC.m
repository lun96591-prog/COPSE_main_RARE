%% ===================================================================
%%  Fig. 3  Paired runs: Q_rw on vs off  (GBC style)
%% -------------------------------------------------------------------
%%  Changes vs the original paired_runs_figure.m:
%%    * sgtitle removed (AGU/GBC: no text labels inside the image)
%%    * Arial lettering; Okabe-Ito colour-blind-safe palette
%%    * panels re-spaced to use the space freed by the title
%%    * 17.8 cm width, 600 dpi, RGB, white background
%%  Data logic and the W_extra integral are unchanged.
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

Q_on  = 0.0608;      % Ea = 42 kJ/mol
Q_off = 0;

w0 = -248.7;  w1 = -237;      % W_extra window
zlo = -265;   zhi = -215;     % zoom range

base = struct();
base.k_rw_si         = 1.0e12;
base.Alk_Si_RW_ratio = 2.66;
base.Si_rw_max       = 15e12;
base.f_si_evo_PTB    = 0.04;

cfg  = [Q_on, Q_off];
runs = struct();
for k = 1:2
    SENS = struct();
    SENS.k_rw_si         = base.k_rw_si;
    SENS.Alk_Si_RW_ratio = base.Alk_Si_RW_ratio;
    SENS.Si_rw_max       = base.Si_rw_max;
    SENS.f_si_evo_PTB    = base.f_si_evo_PTB;
    SENS.Q_rw            = cfg(k);

    fprintf('Running Q_rw = %.4f ...\n', cfg(k));
    r = COPSE_frontend(0); close all

    st = r.state;
    runs(k).Q   = cfg(k);
    runs(k).t   = st.time_myr(:);
    runs(k).T   = st.tempC(:);
    runs(k).DSi = st.DSi_conc_uM(:);
    runs(k).Srw = st.Si_rw(:) / 1e12;
    runs(k).fT  = st.f_T_rw(:);
end

%% ---------------------- derived quantities ------------------------
t_on  = runs(1).t;   t_off = runs(2).t;
T_on  = runs(1).T;   T_off = runs(2).T;

T_off_i = interp1(t_off, T_off, t_on, 'linear', 'extrap');
dT = T_on - T_off_i;

selW    = (t_on >= w0) & (t_on <= w1);
W_extra = trapz(t_on(selW), dT(selW));          % degC . Myr

tw  = t_on(selW);
Tw  = T_on(selW);
[Tpk, ipk] = max(Tw);
t_pk       = tw(ipk);
Tpk_off    = max(T_off_i(selW));

selZ = (t_on >= zlo) & (t_on <= zhi);

idxZ = find(selZ);
[fTmax, jz] = max(runs(1).fT(selZ));
ifT = idxZ(jz);

frac_gt1   = mean(runs(1).fT > 1);
fT_win_min = min(runs(1).fT(selW));
fT_win_max = max(runs(1).fT(selW));

fprintf('\nW_extra (%0.1f to %0.0f Ma) = %.2f degC.Myr\n', w0, w1, W_extra);
fprintf('Window peak: on %.2f degC @ %.2f Ma; off %.2f degC (diff %+.2f)\n', ...
    Tpk, t_pk, Tpk_off, Tpk - Tpk_off);

%% ---------------------- style -------------------------------------
c_on   = [0.000 0.447 0.741];   % Okabe-Ito blue
c_off  = [0.835 0.369 0.000];   % Okabe-Ito vermillion
c_srw  = [0.008 0.620 0.451];   % Okabe-Ito bluish green
c_srw0 = [0.494 0.184 0.556];   % Okabe-Ito purple
c_gry  = [0.800 0.800 0.800];
c_fill = [0.902 0.624 0.000];   % Okabe-Ito orange

xlo = zlo;  xhi = zhi;

figW = 7.0;  figH = 7.6;
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 figW figH]);
set(fig, 'DefaultAxesFontName','Arial', 'DefaultTextFontName','Arial', ...
         'DefaultLegendFontName','Arial', 'DefaultAxesFontSize',8, ...
         'DefaultTextFontSize',8, 'DefaultLegendFontSize',8, ...
         'DefaultAxesLineWidth',0.6, 'DefaultAxesTickDir','in', ...
         'DefaultAxesLayer','top');

ax1 = axes('Parent',fig,'Position',[0.095 0.700 0.845 0.258]); hold(ax1,'on'); box(ax1,'on')
ax2 = axes('Parent',fig,'Position',[0.095 0.410 0.845 0.258]); hold(ax2,'on'); box(ax2,'on')
ax3 = axes('Parent',fig,'Position',[0.095 0.120 0.845 0.258]); hold(ax3,'on'); box(ax3,'on')
ax4 = axes('Parent',fig,'Position',[0.095 0.062 0.845 0.042]); hold(ax4,'on')

%% ---------- (a) GMST with the residual-integration window ----------
y1 = min([T_on; T_off]);  y2 = max([T_on; T_off]);  pad = 0.08*(y2-y1);
xlim(ax1,[xlo xhi]);  ylim(ax1, [y1-pad, y2+pad]);

patch(ax1, [w0 w1 w1 w0], [y1-pad y1-pad y2+pad y2+pad], c_gry, ...
    'FaceAlpha',0.18, 'EdgeColor','none');

patch(ax1, [tw; flipud(tw)], [Tw; flipud(T_off_i(selW))], c_fill, ...
    'FaceAlpha',0.45, 'EdgeColor',[0.80 0.55 0.00], 'LineWidth',0.4);

h1 = plot(ax1, t_on,  T_on,  '-',  'Color', c_on,  'LineWidth', 1.0);
h2 = plot(ax1, t_off, T_off, '--', 'Color', c_off, 'LineWidth', 0.9);
plot(ax1, t_pk, Tpk, 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerSize', 4);
xline(ax1, w0, ':', 'Color',[0.45 0.45 0.45], 'LineWidth',0.6);
xline(ax1, w1, ':', 'Color',[0.45 0.45 0.45], 'LineWidth',0.6);
ylabel(ax1, 'GMST (°C)')
ax1.XTickLabel = {};
text(ax1, t_pk - 1.0, Tpk, sprintf('peak %.2f °C ', Tpk), ...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontSize', 8)
rng1 = (y2+pad) - (y1-pad);
text(ax1, mean([w0 w1]), interp1(t_off, T_off, mean([w0 w1])) - 0.12*rng1, ...
    sprintf('W_{extra} = %.2f °C·Myr', W_extra), ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', ...
    'FontSize', 8, 'Color',[0.60 0.40 0.00], 'FontWeight','bold')
legend(ax1, [h1 h2], {'Q_{rw} > 0', 'Q_{rw} = 0'}, 'Location','southwest', 'FontSize', 8);
grid(ax1,'on')

%% ---------- (b) DSi + Si_rw (dual axis) ---------------------------
yyaxis(ax2,'right')
ax2.YAxis(2).Limits = [0 18];
ax2.YAxis(2).TickValues = 0:3:18;
patch(ax2, [w0 w1 w1 w0], [0 0 18 18], c_gry, 'FaceAlpha',0.18, 'EdgeColor','none');
hb3 = plot(ax2, t_on,  runs(1).Srw, '-.', 'Color', c_srw,  'LineWidth', 1.0);
hb4 = plot(ax2, t_off, runs(2).Srw, ':',  'Color', c_srw0, 'LineWidth', 0.9);
ylabel(ax2, 'Si_{rw} (Tmol Si yr^{-1})')
ax2.YAxis(2).Color = c_srw;

yyaxis(ax2,'left')
hb1 = plot(ax2, t_on,  runs(1).DSi, '-',  'Color', c_on,  'LineWidth', 1.0);
hb2 = plot(ax2, t_off, runs(2).DSi, '--', 'Color', c_off, 'LineWidth', 0.9);
ylabel(ax2, 'DSi (µM)')
ax2.YAxis(1).Color = c_on;

xlim(ax2,[xlo xhi])
xline(ax2, w0, ':', 'Color',[0.45 0.45 0.45], 'LineWidth',0.6);
xline(ax2, w1, ':', 'Color',[0.45 0.45 0.45], 'LineWidth',0.6);
ax2.XTickLabel = {};
legend(ax2, [hb1 hb2 hb3 hb4], ...
    {'DSi, Q_{rw} > 0', 'DSi, Q_{rw} = 0', 'Si_{rw}, Q_{rw} > 0', 'Si_{rw}, Q_{rw} = 0'}, ...
    'Location','northeast', 'FontSize', 8, 'NumColumns', 2);
grid(ax2,'on')

%% ---------- (c) f_T,rw --------------------------------------------
xlim(ax3,[xlo xhi])
yl3 = [min([runs(1).fT(selZ); runs(2).fT(selZ); 1]) - 0.05, ...
       max([runs(1).fT(selZ); runs(2).fT(selZ); 1]) + 0.20];
ylim(ax3, yl3);

patch(ax3, [w0 w1 w1 w0], [yl3(1) yl3(1) yl3(2) yl3(2)], c_gry, ...
    'FaceAlpha',0.18, 'EdgeColor','none');

hf1 = plot(ax3, t_on,  runs(1).fT, '-',  'Color', c_on,  'LineWidth', 1.0);
hf2 = plot(ax3, t_off, runs(2).fT, '--', 'Color', c_off, 'LineWidth', 0.9);
hy  = yline(ax3, 1, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.7);

xline(ax3, w0, ':', 'Color',[0.45 0.45 0.45], 'LineWidth',0.6);
xline(ax3, w1, ':', 'Color',[0.45 0.45 0.45], 'LineWidth',0.6);
ylabel(ax3, 'f_{T,rw} (dimensionless)')
ax3.XTickLabel = {};
legend(ax3, [hf1 hf2 hy], {'Q_{rw} > 0', 'Q_{rw} = 0', 'f_{T,rw} = 1'}, ...
    'Location','northwest', 'FontSize', 8);
grid(ax3,'on')

%% ---------- (d) geological time bar -------------------------------
per = {
    {'Permian'},           -298.9,   -251.902, [0.87 0.84 0.78];
    {'Early','Triassic'},  -251.902, -247.2,   [0.96 0.79 0.69];
    {'Middle','Triassic'}, -247.2,   -237.0,   [0.93 0.86 0.73];
    {'Late','Triassic'},   -237.0,   -201.4,   [0.89 0.89 0.84];
};
xlim(ax4, [xlo xhi]);  ylim(ax4, [0 1]);
for i = 1:size(per,1)
    a0 = max(per{i,2}, xlo);  a1 = min(per{i,3}, xhi);
    if a1 <= a0, continue; end
    patch(ax4, [a0 a1 a1 a0], [0 0 1 1], per{i,4}, ...
        'EdgeColor',[0.32 0.32 0.32], 'LineWidth',0.6);
    text(ax4, 0.5*(a0+a1), 0.5, per{i,1}, 'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', 'FontSize', 8, 'Color',[0.12 0.12 0.12]);
end
ax4.YTick = [];
ylabel(ax4, 'Period', 'FontSize', 8)
xlabel(ax4, 'time (Ma)')

for b = [-251.902, -247.2, -237.0]
    if b > xlo && b < xhi
        xline(ax1, b, ':', 'Color',[0.82 0.82 0.82], 'LineWidth',0.5, 'HandleVisibility','off');
        xline(ax2, b, ':', 'Color',[0.82 0.82 0.82], 'LineWidth',0.5, 'HandleVisibility','off');
        xline(ax3, b, ':', 'Color',[0.82 0.82 0.82], 'LineWidth',0.5, 'HandleVisibility','off');
    end
end

axs = [ax1 ax2 ax3];
lbl = {'(a)','(b)','(c)'};
for k = 1:3
    pp = axs(k).Position;
    annotation(fig, 'textbox', [0.014, pp(2)+pp(4)-0.026, 0.055, 0.026], ...
        'String', lbl{k}, 'FontSize', 9, 'FontWeight', 'bold', ...
        'FontName', 'Arial', 'LineStyle', 'none', ...
        'VerticalAlignment', 'middle');
end

linkaxes([ax1 ax2 ax3 ax4], 'x')
ax4.XTick = xlo:5:xhi;

%% ---------------------- save --------------------------------------
set(findall(fig,'-property','FontName'), 'FontName', 'Arial');
set(fig,'PaperPositionMode','auto','InvertHardcopy','on','Color','w');
print(fig, fullfile(outdir,'fig3_paired_runs_GBC.png'), '-dpng', '-r600');
save( fullfile(outdir,'fig3_paired_runs_GBC.mat'), ...
    'runs','W_extra','t_pk','Tpk','Tpk_off','fTmax','frac_gt1', ...
    'fT_win_min','fT_win_max','w0','w1','Q_on','Q_off','base');

fprintf('\nSaved: %s\n', fullfile(outdir,'fig3_paired_runs_GBC.png'));
