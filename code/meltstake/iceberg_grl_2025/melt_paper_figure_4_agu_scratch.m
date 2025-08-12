% Script to create figure 4 for the melt rate GRL paper. This figure is a
% multi-panel scatter plot comparing melt rate observations to melt rates
% predicted by various parameterizations. Output for parameterizations is
% created using the script "ms_melt_models.m"
%
% KJW
% 7 Apr 2024
set(0,'defaulttextinterpreter','latex');

clear

% load data
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

% sections
sections = [8 3 29 9];
lbls = {'A','B','C','D'};
mkrs = {'>','^','square','hexagram'};
lgd_lbls = {'A (crossflow)','B (plume)','C (wave)','D (eddy)'};

% observations
m_obs = msTable2Vector(ms_tbl.m);
m_ci = ms_tbl.m_ci;
m_obs = m_obs*.24; % cm/hr --> m/day
m_ci = m_ci*.24;
[u,u_ci] = msTable2Vector(ms_tbl.u0);
up = ms_tbl.uprimesqrd;
[T,T_ci] = msTable2Vector(ms_tbl.T);

m_obs(end) = nan;
m_ci(end) = nan;

% models
load ../../data/melt_models.mat

% melt error
r = repmat(m_obs,[1 size(m_mean,2)])./m_mean;
r_ci = repmat(m_ci,[1 size(m_mean,2)])./m_mean;

% UT contours
n_cont = 100;
u_cont = linspace(0,.15,n_cont);
t_cont = linspace(0,11,n_cont);
[U_cont,T_cont] = meshgrid(u_cont,t_cont);
UT_cont = U_cont.*T_cont;
m_cont = solve3Eqn(U_cont,T_cont,25,0)*86400;
cont_lvls = .05:.15:1.5;

%% plot
fs  = 11;
lw = 1;
clear ax
clear h

% padding
pxi = .11;
pxo = .08;
pyi = .05;
pyo = .16;
pad = [.11 .08 .05 .16];
shift = [.02 0];

%%% 2 panel figure %%%%
% create and size figure
figsize = [18 9.5];
fig1 = figure(4); clf
setFigureSize(fig1,figsize);

model_lbls = {'3-equation','convective','modified 3-eqn','bulk iceberg','3-eqn (coarse)'};
title_lbls = {'','','(e)','(d)','(c)'};
model_num = 1;

% create axes
ax1 = axes(fig1,'position',axgridpos(1,2,1,pad,shift) - [0 .04 0 0]);
hold on
box on

% plot
% panel 1: melt error vs U and T
pt_lbls = upper({'a','b','c','d'});
marker_scale = 75;
contour(ax1,u_cont,t_cont,m_cont,cont_lvls,'k:')
% errorbar(ax1,u,T,-T_ci,T_ci,-u_ci,u_ci,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
scatter(ax1,u,T,marker_scale,r(:,model_num),'o','filled','markeredgecolor','k')
text(ax1,.15,.98,[model_lbls{model_num}],'units','normalized','fontsize',fs,...
        'verticalalignment','top','horizontalalignment','left')
for i = 1:length(sections)
    scatter(ax1,u(sections(i)),T(sections(i)),marker_scale,r(sections(i),model_num),'o','filled','markeredgecolor','r','linewidth',lw)
    text(ax1,u(sections(i))-.005,T(sections(i))+.005,sprintf('%s',pt_lbls{i}),'color','k',...
        'horizontalalignment','center','verticalalignment','bottom',...
        'fontsize',fs+2,'backgroundcolor','none')
end

axpos = ax1.Position;
cbar = colorbar(ax1);
cbar.Location = 'northoutside';
box on

xlim(ax1,[0 .13])
ylim(ax1,[2.5 11])
clim(ax1,[0 3.5])

% labels
xlabel(ax1,'$U$ [m/s]','fontsize',fs)
ylabel(ax1,'$T$ [$^\circ$C]','fontsize',fs)
cbar.Label.String = '$m_\mathrm{obs}/m_\mathrm{mod}$';
cbar.Label.FontSize = fs;
cbar.Label.Interpreter = 'latex';
cmocean('balance','pivot',1)

% tighten up colorbar spacing
ax1.Position = axpos;
cbar.Position(4) = .04;
cbar.Position(2) = sum(ax1.Position([2 4])) + .012;

% panel 2: obs vs model
marker_clr = 0.5*[1 1 1];
ax2 = axes(fig1,'position',axgridpos(1,2,2,pad,shift) - [0 .04 0 0]);
hold on
box on
plot(ax2,[0 15],1*[0 15],'k-')
errorbar(ax2,m_mean(:,model_num),m_obs,-m_ci,m_ci,...
    'color',.5*[1 1 1],'linestyle','none','linewidth',lw)
plot(ax2,m_mean(:,model_num),m_obs,'ko','markersize',5,'markerfacecolor',marker_clr)
xlim([0 2])
ylim([0 2])

ylabel('observed melt rate [m/day]','fontsize',fs)
xlabel('modeled melt rate [m/day]','fontsize',fs)


text(ax2,.01,.98,['(b) ' model_lbls{model_num}],'units','normalized','fontsize',fs,...
    'verticalalignment','top','horizontalalignment','left')

% add sections
for j = 1:length(sections)
    h(j) = plot(ax2,m_mean(sections(j),model_num),m_obs(sections(j)),'ko','markersize',7,...
        'markerfacecolor',colors(j),'marker',mkrs{j});
end
if model_num == 1
    lgd = legend(ax2,h,lgd_lbls,'location','southeast','fontsize',fs-1,'color','w');
    lgd.Position = lgd.Position + [.015 -.035 0 0];
    lgd.Interpreter = 'latex';
end

if model_num==3
    xlim([0 13])
    ylim([0 13])
%     lgd.Location = 'best';
end

% histogram inset
axes(fig1,'position',[0.85 0.55 0.11 0.2],'units','normalized');
hold on
box on
histogram(log10(r(:,1)),-1.5:0.15:2.3,'edgecolor','none','facecolor',0.3*[1 1 1])
histogram(log10(r(:,model_num)),-1.5:0.15:2.3,'facecolor',colors(7),'facealpha',1)
plot(log10(median(r(:,model_num),'omitnan'))*[1 1],[0 15],'k--','linewidth',1)
xlim([-2 2])
ylim([0 15])
set(gca,'xtick',-2:2)
set(gca,'xticklabelrotation',0)

xlabel('$\log_{10}\big(\frac{m_{\mathrm{obs}}}{m_{\mathrm{mod}}}\big)$','fontsize',fs+1)
ylabel('occur.','fontsize',fs)

%%% 1 panel figure %%%%
% create and size figure
figsize = [6 6];
fig2 = figure(5); clf
setFigureSize(fig2,figsize);

% panel 2: obs vs model
ax3 = axes(fig2);
hold on
box on
plot(ax3,[0 15],1*[0 15],'k-')
errorbar(ax3,m_mean(:,model_num),m_obs,-m_ci,m_ci,...
    'color',.5*[1 1 1],'linestyle','none','linewidth',lw)
plot(ax3,m_mean(:,model_num),m_obs,'ko','markersize',5,'markerfacecolor',marker_clr)
xlim([0 2])
ylim([0 2])

ylabel('observed melt rate [m/day]','fontsize',fs)
xlabel('modeled melt rate [m/day]','fontsize',fs)

% title([title_lbls{model_num} ' ' model_lbls{model_num}],'fontsize',fs)
% text(ax3,.01,.98,[title_lbls{model_num} ' ' model_lbls{model_num}],'units','normalized','fontsize',fs,...
%     'verticalalignment','top','horizontalalignment','left')
text(ax3,.15,.98,[model_lbls{model_num}],'units','normalized','fontsize',fs,...
    'verticalalignment','top','horizontalalignment','left')

% add sections
for j = 1:length(sections)
    h(j) = plot(ax3,m_mean(sections(j),model_num),m_obs(sections(j)),'ko','markersize',7,...
        'markerfacecolor',colors(j),'marker',mkrs{j});
end
if model_num == 1
%     lgd = legend(ax3,h,lgd_lbls,'location','southeast','fontsize',fs,'color','w');
end

if model_num==3
    xlim([0 13])
    ylim([0 13])
%     lgd.Location = 'best';
end

% histogram inset
ax4 = axes(fig2,'position',[0.56 0.49 0.3 0.32],'units','normalized');
hold on
box on
histogram(log10(r(:,1)),-1.5:0.15:2.3,'edgecolor','none','facecolor',0.3*[1 1 1])
histogram(log10(r(:,model_num)),-1.5:0.15:2.3,'facecolor',colors(7),'facealpha',1)
plot(log10(median(r(:,model_num),'omitnan'))*[1 1],[0 15],'k--','linewidth',1)
xlim([-2 2])
ylim([0 15])
set(ax4,'clipping','off')
set(ax4,'xtick',-2:2)
set(ax4,'xticklabelrotation',0)

xlabel('$\log_{10}\big(\frac{m_{\mathrm{obs}}}{m_{\mathrm{mod}}}\big)$','fontsize',fs+1)
ylabel('occur.','fontsize',fs)

% % create and size figure
% figsize = [6 6];
% fig2 = figure(5); clf
% setFigureSize(fig2,figsize);
% 
% % obs vs model
% ax3 = axes(fig2);
% hold on
% box on
% plot(ax3,[0 15],1*[0 15],'k-')
% errorbar(ax3,m_mean(:,model_num),m_obs,-m_ci,m_ci,...
%     'color',.5*[1 1 1],'linestyle','none','linewidth',lw)
% plot(ax3,m_mean(:,model_num),m_obs,'ko','markersize',5,'markerfacecolor',marker_clr)
% xlim([0 2])
% ylim([0 2])
% 
% ylabel('observed melt rate [m/day]','fontsize',fs)
% xlabel('modeled melt rate [m/day]','fontsize',fs)
% 
% % title([title_lbls{model_num} ' ' model_lbls{model_num}],'fontsize',fs)
% text(ax3,.01,.98,model_lbls{model_num},'units','normalized','fontsize',fs,...
%     'verticalalignment','top','horizontalalignment','left')
% 
% if model_num==3
%     xlim([0 13])
%     ylim([0 13])
% %     lgd.Location = 'best';
% end

%%% histogram %%%%
% % create and size figure
% figsize = [3 3];
% fig3 = figure(6); clf
% setFigureSize(fig3,figsize);
% 
% hold on
% box on
% % hist_mods = [1 3 4 5];
% % for i = hist_mods
% %     histogram(log10(r(:,i)),-1.5:0.1:2.3)
% % end
% histogram(log10(r(:,1)),-1.5:0.15:2.3,'edgecolor','none','facecolor',0.3*[1 1 1])
% histogram(log10(r(:,model_num)),-1.5:0.15:2.3,'facecolor',colors(7),'facealpha',1)
% plot(log10(median(r(:,model_num),'omitnan'))*[1 1],[0 15],'k--','linewidth',1)
% xlim([-2 2])
% ylim([0 15])
% set(gca,'clipping','off')
% set(gca,'xtick',-2:2)
% 
% xlabel('$\log_{10}\big(\frac{m_{obs}}{m_{mod}}\big)$','fontsize',fs)
% ylabel('occur.','fontsize',fs)
% % hist_lgd = legend(model_lbls(hist_mods));