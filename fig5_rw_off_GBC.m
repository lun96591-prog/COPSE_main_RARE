%% ===================================================================
%%  Fig. 5  Reverse weathering as a process control  (GBC style)
%% -------------------------------------------------------------------
%%  Changes vs the original fig4_rw_off.m:
%%    * sgtitle removed (AGU/GBC: no labels inside the image)
%%    * Arial lettering; Okabe-Ito palette (blue / vermillion /
%%      bluish green) replaces the red-green pairing, which is unsafe
%%      for deuteranope readers
%%    * legend kept but re-anchored in the space freed by the title
%%    * 17.8 cm width, 600 dpi, RGB, white background
%%  Data logic and the three scenarios are unchanged.
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

base = struct('k_rw_si', 1.0e12, 'Q_rw', 0.0608, ...
              'Alk_Si_RW_ratio', 2.66, 'Si_rw_max', 15e12);

scen = {
    'calibrated',                base.k_rw_si, base.Q_rw;
    'no reverse weathering',     0.0,           base.Q_rw;
    'no temperature dependence', base.k_rw_si,  0.0;
};
nS = size(scen,1);

run_data = struct('t',{},'pC',{},'T',{},'TAlk',{});

for s = 1:nS
    fprintf('Running scenario %d/%d: %s ...\n', s, nS, scen{s,1})
    clear global SENS Gtune Ctune PYRtune GYPtune Atune Otune Stune
    global SENS Gtune
    Gtune = [];

    SENS = struct();
    SENS.k_rw_si         = scen{s,2};
    SENS.Q_rw            = scen{s,3};
    SENS.Alk_Si_RW_ratio = base.Alk_Si_RW_ratio;
    SENS.Si_rw_max       = base.Si_rw_max;

    rr = COPSE_frontend(0); close all
    st = rr.state;

    run_data(s).t    = st.time_myr(:);
    run_data(s).pC   = st.RCO2(:)*280;
    run_data(s).T    = st.tempC(:);
    run_data(s).TAlk = st.Alk_conc_mM(:);
end

for s = 1:nS
    t   = run_data(s).t;
    iET = (t >= -253) & (t <= -248);
    iMo = (t >= -1)   & (t <= 0);
    fprintf('%-26s ET: pCO2=%6.0f ppm  T=%5.2f °C  TAlk=%.2f | modern: pCO2=%5.0f ppm\n', ...
        scen{s,1}, mean(run_data(s).pC(iET)), mean(run_data(s).T(iET)), ...
        mean(run_data(s).TAlk(iET)), mean(run_data(s).pC(iMo)));
end

%% ---------------------- style -------------------------------------
cAll = [0.000 0.447 0.741;     % A calibrated        (Okabe-Ito blue)
        0.835 0.369 0.000;     % B rw off            (Okabe-Ito vermillion)
        0.008 0.620 0.451];    % C no T dependence   (Okabe-Ito bluish green)
sty  = {'-','--','-.'};
lw   = [1.0 0.9 0.9];
c_gry = [0.80 0.80 0.80];

xlo = -260;  xhi = -237;

figW = 7.0;  figH = 4.85;
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 figW figH]);
set(fig, 'DefaultAxesFontName','Arial', 'DefaultTextFontName','Arial', ...
         'DefaultLegendFontName','Arial', 'DefaultAxesFontSize',8, ...
         'DefaultTextFontSize',8, 'DefaultLegendFontSize',8, ...
         'DefaultAxesLineWidth',0.6, 'DefaultAxesTickDir','in', ...
         'DefaultAxesLayer','top');

ax1 = axes('Parent',fig,'Position',[0.098 0.660 0.855 0.215]); hold(ax1,'on'); box(ax1,'on')
ax2 = axes('Parent',fig,'Position',[0.098 0.420 0.855 0.190]); hold(ax2,'on'); box(ax2,'on')
ax3 = axes('Parent',fig,'Position',[0.098 0.180 0.855 0.190]); hold(ax3,'on'); box(ax3,'on')
ax4 = axes('Parent',fig,'Position',[0.098 0.086 0.855 0.046]); hold(ax4,'on');

axs  = [ax1 ax2 ax3];
vars = {'pC','T','TAlk'};
labs = {'pCO_2 (ppm)', 'GMST (°C)', 'TAlk (mmol kg^{-1})'};

h = gobjects(nS,3);
for k = 1:3
    a = axs(k);
    for s = 1:nS
        h(s,k) = plot(a, run_data(s).t, run_data(s).(vars{k}), sty{s}, ...
            'Color', cAll(s,:), 'LineWidth', lw(s));
    end
    xlim(a, [xlo xhi])
    yl = ylim(a);
    ylim(a, [max(0, yl(1)) yl(2)*1.08])
    y2 = ylim(a);
    hp = patch(a, [-253 -248 -248 -253], [y2(1) y2(1) y2(2) y2(2)], c_gry, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    uistack(hp, 'bottom')
    ylim(a, y2)

    ylabel(a, labs{k})
    grid(a, 'on')
    a.XTick = xlo:5:xhi;
    a.XTickLabel = {};
end

lg = legend(ax1, [h(1,1) h(2,1) h(3,1)], ...
    {'calibrated', 'no reverse weathering (k_{rw} = 0)', ...
     'no temperature dependence (Q_{rw} = 0)'}, ...
    'Orientation','horizontal', 'Location','northoutside', 'FontSize', 8);
lg.Box = 'off';

%% ---------- geological time bar -----------------------------------
per = {
    {'Permian'},           -298.9,   -251.902, [0.87 0.84 0.78];
    {'Early','Triassic'},  -251.902, -247.2,   [0.96 0.79 0.69];
    {'Middle','Triassic'}, -247.2,   -237.0,   [0.93 0.86 0.73];
    {'Late','Triassic'},   -237.0,   -201.4,   [0.89 0.89 0.84];
};
xlim(ax4, [xlo xhi]);  ylim(ax4, [0 1]);  box(ax4, 'on')
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
ax4.XTick = xlo:5:xhi;
xlabel(ax4, 'time (Ma)')

for b = [-251.902, -247.2]
    if b > xlo && b < xhi
        xline(ax1, b, ':', 'Color',[0.82 0.82 0.82], 'LineWidth',0.5, 'HandleVisibility','off');
        xline(ax2, b, ':', 'Color',[0.82 0.82 0.82], 'LineWidth',0.5, 'HandleVisibility','off');
        xline(ax3, b, ':', 'Color',[0.82 0.82 0.82], 'LineWidth',0.5, 'HandleVisibility','off');
    end
end

lbl = {'(a)','(b)','(c)'};
for k = 1:3
    pp = axs(k).Position;
    annotation(fig, 'textbox', [0.010, pp(2)+pp(4)-0.020, 0.040, 0.024], ...
        'String', lbl{k}, 'FontSize', 9, 'FontWeight','bold', ...
        'FontName','Arial', 'LineStyle','none', ...
        'VerticalAlignment','middle');
end

linkaxes([ax1 ax2 ax3 ax4], 'x')
ax4.XTick = xlo:5:xhi;

%% ---------------------- save --------------------------------------
set(findall(fig,'-property','FontName'), 'FontName', 'Arial');
set(fig,'PaperPositionMode','auto','InvertHardcopy','on','Color','w');
print(fig, fullfile(outdir,'fig5_rw_off_GBC.png'), '-dpng', '-r600');
save(fullfile(outdir,'fig5_rw_off_GBC.mat'), 'scen','base','run_data','xlo','xhi');

fprintf('\nSaved: %s\n', fullfile(outdir,'fig5_rw_off_GBC.png'));
