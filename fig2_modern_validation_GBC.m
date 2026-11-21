%% ===================================================================
%%  Fig. 2  Modern-state calibration of the composite model  (GBC style)
%% -------------------------------------------------------------------
%%  Changes vs the original fig5_modern_validation.m:
%%    * all in-figure titles removed (AGU/GBC: no labels or captions
%%      inside the image file)
%%    * Arial lettering (AGU recommends a sans-serif face)
%%    * Okabe-Ito colour-blind-safe palette
%%    * 17.8 cm double-column width, 600 dpi, RGB, white background
%%    * panel labels (a)/(b) moved outside the axes, 9 pt bold
%%  Data logic is unchanged.
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

SENS = struct();
SENS.k_rw_si         = 1.0e12;
SENS.Alk_Si_RW_ratio = 2.66;
SENS.Si_rw_max       = 15e12;
SENS.f_si_evo_PTB    = 0.04;
SENS.Q_rw            = 0.0608;

rr = COPSE_frontend(0); close all
st = rr.state;
t  = st.time_myr(:);
[~, j0] = min(abs(t));          % instantaneous modern state at t = 0

mod.pCO2 = st.RCO2(j0)*280;
mod.DSi  = st.DSi_conc_uM(j0);
mod.TAlk = st.Alk_conc_mM(j0);
mod.pH   = st.pH(j0);
mod.T    = st.tempC(j0);
mod.om   = st.Omega_cal(j0);
mod.silw = st.Si_silw(j0)/1e12;
mod.sfw  = st.Si_sfw(j0)/1e12;
mod.bio  = st.Si_bio(j0)/1e12;
mod.rw   = st.Si_rw(j0)/1e12;

fprintf('Modern state: pCO2=%.1f ppm | DSi=%.1f | TAlk=%.2f | pH=%.3f | T=%.2f\n', ...
    mod.pCO2, mod.DSi, mod.TAlk, mod.pH, mod.T);
fprintf('Si fluxes: silw=%.2f | sfw=%.2f | bio=%.2f | rw=%.2f Tmol/yr\n', ...
    mod.silw, mod.sfw, mod.bio, mod.rw);

%% ---------------------- style -------------------------------------
cB    = [0.000 0.447 0.741];   % Okabe-Ito blue
cV    = [0.835 0.369 0.000];   % Okabe-Ito vermillion
cBand = [0.850 0.900 0.850];   % +/-20 % band
cGry  = [0.550 0.550 0.550];

figW = 7.0;  figH = 3.2;
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 figW figH]);
set(fig, 'DefaultAxesFontName','Arial', 'DefaultTextFontName','Arial', ...
         'DefaultLegendFontName','Arial', 'DefaultAxesFontSize',8, ...
         'DefaultTextFontSize',8, 'DefaultLegendFontSize',8, ...
         'DefaultAxesLineWidth',0.6, 'DefaultAxesTickDir','in', ...
         'DefaultAxesLayer','top');

%% ---------------- (a) state variables: model / observation -------
axA = axes('Parent',fig,'Position',[0.110 0.195 0.390 0.740]); hold(axA,'on'); box(axA,'on')

lab = {'pCO_{2}', 'DSi', 'TAlk', 'pH', 'GMST'};
mo  = [mod.pCO2, mod.DSi, mod.TAlk, mod.pH, mod.T];
ob  = [280, 80, 2.30, 8.05, 15.0];
rat = mo ./ ob;
n   = numel(lab);

hp = patch(axA, [0.8 1.2 1.2 0.8], [0 0 n+0.6 n+0.6], cBand, ...
    'FaceAlpha',0.55, 'EdgeColor','none');
uistack(hp,'bottom')

for i = 1:n
    y = n + 1 - i;
    plot(axA, [1 rat(i)], [y y], '-', 'Color', cGry, 'LineWidth', 1.0);
    plot(axA, rat(i), y, 'o', 'MarkerSize', 4.5, ...
        'MarkerFaceColor', cB, 'MarkerEdgeColor', cB);
    if rat(i) < 1
        text(axA, rat(i)*0.955, y, sprintf('%.2f', rat(i)), ...
            'FontSize', 7.5, 'Color', cB, 'VerticalAlignment','middle', ...
            'HorizontalAlignment','right');
    else
        text(axA, rat(i)*1.045, y, sprintf('%.2f', rat(i)), ...
            'FontSize', 7.5, 'Color', cB, 'VerticalAlignment','middle', ...
            'HorizontalAlignment','left');
    end
end
xline(axA, 1, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);
set(axA, 'XScale','log', 'XLim',[0.8 1.8], 'YLim',[0.5 n+0.6], ...
    'YTick', 1:n, 'YTickLabel', flip(lab))
xlabel(axA, 'model / observation')
grid(axA, 'on')

%% ---------------- (b) modern silicon fluxes: model / published ----
axB = axes('Parent',fig,'Position',[0.620 0.195 0.360 0.740]); hold(axB,'on'); box(axB,'on')

fl   = {'F_{silw}', 'F_{bio}', 'Si_{rw}'};
fmod = [mod.silw, mod.bio, mod.rw];
fpub = [7.23, 5.58, 2.00];      % Treguer et al. (2021); Rahman et al. (2017)
frat = fmod ./ fpub;
nf   = numel(fl);

hp = patch(axB, [0.8 1.2 1.2 0.8], [0 0 nf+0.6 nf+0.6], cBand, ...
    'FaceAlpha',0.55, 'EdgeColor','none');
uistack(hp,'bottom')

for i = 1:nf
    y = nf + 1 - i;
    plot(axB, [1 frat(i)], [y y], '-', 'Color', cGry, 'LineWidth', 1.0);
    plot(axB, frat(i), y, 'o', 'MarkerSize', 4.5, ...
        'MarkerFaceColor', cV, 'MarkerEdgeColor', cV);
    if frat(i) < 1
        text(axB, frat(i)*0.955, y, sprintf('%.2f', frat(i)), ...
            'FontSize', 7.5, 'Color', cV, 'VerticalAlignment','middle', ...
            'HorizontalAlignment','right');
    else
        text(axB, frat(i)*1.045, y, sprintf('%.2f', frat(i)), ...
            'FontSize', 7.5, 'Color', cV, 'VerticalAlignment','middle', ...
            'HorizontalAlignment','left');
    end
end
xline(axB, 1, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);
set(axB, 'XScale','log', 'XLim',[0.8 1.8], 'YLim',[0.5 nf+0.6], ...
    'YTick', 1:nf, 'YTickLabel', flip(fl))
xlabel(axB, 'model / published')
grid(axB, 'on')

%% ---------------------- panel labels ------------------------------
posA = axA.Position;  posB = axB.Position;
annotation(fig,'textbox',[posA(1)-0.048, posA(2)+posA(4)+0.012, 0.05, 0.03], ...
    'String','(a)','FontSize',9,'FontWeight','bold','FontName','Arial', ...
    'LineStyle','none','VerticalAlignment','middle');
annotation(fig,'textbox',[posB(1)-0.048, posB(2)+posB(4)+0.012, 0.05, 0.03], ...
    'String','(b)','FontSize',9,'FontWeight','bold','FontName','Arial', ...
    'LineStyle','none','VerticalAlignment','middle');

%% ---------------------- save --------------------------------------
set(findall(fig,'-property','FontName'), 'FontName', 'Arial');
set(fig,'PaperPositionMode','auto','InvertHardcopy','on','Color','w');
print(fig, fullfile(outdir,'fig2_modern_validation_GBC.png'), '-dpng', '-r600');
save(fullfile(outdir,'fig2_modern_validation_GBC.mat'), 'mod','rat','frat','lab','fl');

fprintf('\nSaved: %s\n', fullfile(outdir,'fig2_modern_validation_GBC.png'));
