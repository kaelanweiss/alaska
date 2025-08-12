% Script to create figure 2 for the melt rate GRL paper. This figure has a
% panel showing the full observational parameter space as a scatter plot
% with dimensions (u, T, m, S) where the last two correspond to size and
% color. An additional 3 panels show melt rate as a function of u, T, and S
% separately.
%
% KJW
% 7 Apr 2024
set(0,'defaulttextinterpreter','latex');

clear

addpath('..')

% load data
ms_tbl = loadMSInfo(1:19,'manualwindows');

% cm/hr to m/dy
cmh2md = 0.24;

% extract fields
[m,~] = msTable2Vector(ms_tbl.m);
[T,T_ci] = msTable2Vector(ms_tbl.T);
[S,S_ci] = msTable2Vector(ms_tbl.S);
[u,u_ci] = msTable2Vector(ms_tbl.u_max);
m_ci = ms_tbl.m_ci;

% convert melt rate to m/day
m = m*cmh2md;
m_ci = m_ci*cmh2md;
m(end) = nan;

% standard error
tau_decorr = ms_tbl.tau_decorr;
Neff = 60*ms_tbl.Duration./tau_decorr;
NS = ms_tbl.nS;

T_std = T_ci;
u_std = u_ci;
S_std = S_ci;

T_ci = T_ci./sqrt(Neff);
S_ci = S_ci./sqrt(NS);
u_ci = u_ci./sqrt(Neff);

%% UT contours
n_cont = 100;
u_cont = linspace(0,.15,n_cont);
t_cont = linspace(0,11,n_cont);
[U_cont,T_cont] = meshgrid(u_cont,t_cont);
UT_cont = U_cont.*T_cont;
m_cont = solve3Eqn(U_cont,T_cont,25,0)*86400;
cont_lvls = .05:.15:1.5;

%% plot
lw = 1;
fs = 11;
m_ex = [.2 .5 1 2];
marker_scale = 75;
ldg_lbls = cell(1,length(m_ex));
dy_lgd = 0.4;

pts = [8 3 29 9];
pt_lbls = upper({'a','b','c','d'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 2.1 parameter space (u vs T with colors corresponding to S)
figure(21); clf
% create and size figure
figsize = [10 10];
setFigureSize(figure(21),figsize);

% plot
ax1 = axes(figure(21),'fontsize',fs-2);
hold on
contour(u_cont,t_cont,m_cont,cont_lvls,'k:')
errorbar(u,T,-T_ci,T_ci,-u_ci,u_ci,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
scatter(u,T,marker_scale,m,'o','filled','markeredgecolor','k')
text(.01,.98,'(a)','units','normalized','fontsize',fs,...
        'verticalalignment','top','horizontalalignment','left','interpreter','tex')
for i = 1:length(pts)
    scatter(u(pts(i)),T(pts(i)),marker_scale,m(pts(i)),'o','filled','markeredgecolor','r','linewidth',lw)
    text(ax1,u(pts(i))-.005,T(pts(i))+.005,sprintf('%s',pt_lbls{i}),'color','k',...
        'horizontalalignment','center','verticalalignment','bottom',...
        'fontsize',fs+2,'backgroundcolor','none')
end

cbar = colorbar;
cbar.Location = 'northoutside';
box on

xlim([0 .13])
ylim([2.5 11])
clim([.2 1.5])

% labels
xlabel('$\overline{\textsf{U}_\textsf{max}}$ \textsf{[m/s]}','fontsize',fs)
ylabel('$\overline{\textsf{T}}$ \textsf{[K]}','fontsize',fs)
cbar.Label.String = 'm_{obs} [m/day]';
cbar.Label.FontSize = fs;
cmocean('-ice')

% contour text label
text(0.108,7.2,'$\overline{\textsf{U}_\textsf{max}}\cdot\overline{\textsf{T}}\textsf{=}$\textsf{const.}','rotation',-45,...
    'verticalalignment','middle','horizontalalignment','center')

% tighten up colorbar spacing
axpos = ax1.Position;
ax1.Position = axpos + [0 0 0 .05];
cbar.Position(4) = .04;
cbar.Position(2) = sum(ax1.Position([2 4])) + .012;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 2.2 (m vs u,T,S)
figure(22); clf
% create and size figure
figsize = [9 11];
setFigureSize(figure(22),figsize);

x_data = {u,T,S};
x_ci = {u_ci,T_ci,S_ci};
x_lbls = {'$\overline{\textsf{U}_\textsf{max}}$ \textsf{[m/s]}','$\overline{\textsf{T}}$ \textsf{[K]}','$\overline{\textsf{S}_\textsf{w}}$ \textsf{[psu]}'};
c_data = [2,1,2];
c_cmap = {'-solar','tempo','-solar'};
c_lbls = {'$\overline{\textsf{T}}$ \textsf{[K]}','$\overline{\textsf{U}_\textsf{max}}$ \textsf{[m/s]}','$\overline{\textsf{T}}$ \textsf{[K]}'};
c_lims = {[0 .12],[3 10]};
clrs = [colors(1); colors(3); colors(2)];
% pnl_lbls = {'$m_\mathrm{obs}$ vs. $\overline{U_\mathrm{max}}$','$m_\mathrm{obs}$ vs. $\overline{T}$','$m_\mathrm{obs}$ vs. $\overline{S}$'};
for i = 1:3
    axes(figure(22),'position',axgridpos(3,1,i,[.1,.12,.15,.1],[-.04 0]),'fontsize',fs-2)
    hold on
    errorbar(x_data{i},m,-m_ci,m_ci,-x_ci{i},x_ci{i},'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
    scatter(x_data{i},m,18,x_data{c_data(i)},'filled','o','markeredgecolor','k')
    xlabel(x_lbls{i},'fontsize',fs)
%     ylabel('$\textsf{m}_\textsf{obs}$ \textsf{[m/day]}','fontsize',fs)
    ylabel('m_{obs} [m/day]','fontsize',fs,'interpreter','tex')
    ylim([0 2.2])
    box on
    clim(c_lims{c_data(i)})
    if i==3
        xlim([22 29])
        set(gca,'clipping','off')
    end
    text(.01,.98,['(' char(97+i) ')'],'units','normalized','fontsize',fs,...
        'verticalalignment','top','horizontalalignment','left','interpreter','tex')
    cbar = colorbar;
    cbar.Position = cbarpos(gca,.005,.02);
    cbar.Label.String = c_lbls{i};
    cbar.Label.FontSize = fs;
    cbar.Label.Interpreter = 'latex';
    cmocean(c_cmap{i})
end