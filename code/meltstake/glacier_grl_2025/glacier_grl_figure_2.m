% Script to create panels for figure 2 of the glacier melt rate paper
% (chapter 2). This figure contains 3 scatter plots, all with m_obs on the
% y axis. The x axes of the 3 panels are velocity scale, m_3eqn, and
% m_3eqn,coarse.
%
% KJW
% 22 Sep 2025

clear

addpath('..')
load glacier_clrs

% load data and keep track of iceberg/glacier deps
ms_tbl = loadMSInfo('segments');
inc = ms_tbl.include;
idx_berg = ms_tbl.Number <= 19 & inc;
idx_glac = false(size(ms_tbl,1),3);
for i = 1:3
    idx_glac(:,i) = ms_tbl.Number == (i+25) & inc;
end

% load model results
mod_glac = load('../../../data/melt_models_glacier.mat');
mod_berg = load('../../../data/melt_models.mat');
mod_glac.m_mean = mod_glac.m_mean(any(idx_glac(ms_tbl.Start>datetime(2024,7,1),:),2),:);
mod_berg.m_mean = mod_berg.m_mean(idx_berg(ms_tbl.Start<datetime(2024,7,1)),:);

% cm/hr to m/dy
cmh2md = 0.24;

% extract fields
[m,~] = msTable2Vector(ms_tbl.m);
[T,T_ci] = msTable2Vector(ms_tbl.T);
[u,u_ci] = msTable2Vector(ms_tbl.u_max);
[m_ci_l,m_ci_h] = msTable2Vector(ms_tbl.m_ci);

% convert melt rate to m/day
m = m*cmh2md;
m_ci_l = m_ci_l*cmh2md;
m_ci_h = m_ci_h*cmh2md;

% melt ratio
r = m(inc)./[mod_berg.m_mean(:,1); mod_glac.m_mean(:,1)];
r_crs = m(inc)./[mod_berg.m_mean(:,5); mod_glac.m_mean(:,5)];

% standard error
tau_decorr = ms_tbl.tau_decorr;
Neff = 60*ms_tbl.Duration./tau_decorr;

T_ci = T_ci./sqrt(Neff);
u_ci = u_ci./sqrt(Neff);

%% plot
lw = 0.5;
fs = 11;
marker_scale = 40;
berg_marker = 'diamond';

hist_lvls = -1.5:0.15:2.3;
hist_loc = [.172 .04];
hist_size = [0.09 0.24];

pad = [.035 .0 .06 .18];
shift = [0.35*pad(3) .04];

fig = figure(1); clf
setFigureSize(fig,[15 6])

clear ax
for i = 1:3
    ax(i) = axes(fig,'position',axgridpos(1,3,i,pad,shift),'fontsize',fs-2);
    box on
    hold on
end

dep_lbls = cell(1,3);
for i = 1:3
    dep_lbls{i} = sprintf('%dm',round(mean(ms_tbl.depth(idx_glac(:,i)))));
end

% panel 1: velocity and temperature
clear h1
errorbar(ax(1),u,m,m_ci_l,m_ci_h,-u_ci,u_ci,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
scatter(ax(1),u(idx_berg),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
for i = 1:size(idx_glac,2)
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    h1(i) = scatter(ax(1),u(idx_glac(:,i)),m(idx_glac(:,i)),ms,T(idx_glac(:,i)),smbs{i},'filled','markeredgecolor','k');
end
cmocean('sol',ax(1))
% ax(1).Position(3) = 0.92*ax(1).Position(3);
cbar = colorbar(ax(1),'orientation','horizontal','location','manual');
cbar.Position = [0.22 0.28 0.11 0.04];
cbar.Position = [ax(1).Position(1) sum(ax(1).Position([2 4]))+.01 ax(1).Position(3) .05];
% cbar.Label.String = 'T [K]';
clim(ax(1),[5 6.5])
text(ax(1),cbar.Position(1)-.25,cbar.Position(2)+.1,{'(T-T_f)','[K]'},'fontsize',fs-2,'horizontalalignment','center','verticalalignment','bottom','units','normalized')

% panel 2: m_3eqn
plot(ax(2),[0 10],[0 10],'k-','linewidth',lw)
plot(ax(2),[0 10],1.5*[0 10],'k--','linewidth',lw)
errorbar(ax(2),[mod_berg.m_mean(:,1); mod_glac.m_mean(:,1)],m(inc),m_ci_l(inc),m_ci_h(inc),'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
scatter(ax(2),mod_berg.m_mean(:,1),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
for i = 1:size(idx_glac,2)
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    scatter(ax(2),mod_glac.m_mean(idx_glac(find(any(idx_glac,2),1):end,i),1),m(idx_glac(:,i)),ms,dep_clrs(i,:),smbs{i},'filled','markeredgecolor','k');
end
% histogram
ax_in2 = axes(fig,'position',[ax(2).Position(1:2)+hist_loc 0 0]+[0 0 hist_size],'units','normalized','fontsize',fs-3,'xaxislocation','top');
hold(ax_in2,'on')
box(ax_in2,'on')
histogram(ax_in2,log10(r),hist_lvls,'facecolor',0.*[1 1 1],'facealpha',1)
histogram(ax_in2,log10(r(any(idx_glac(inc,:),2))),hist_lvls,'facecolor',0.7*[1 1 1],'facealpha',1)
% plot(ax_in2,log10(median(r,'omitnan'))*[1 1],[0 15],'k-','linewidth',1)
xlim(ax_in2,[-0.5 1.1])
ylim(ax_in2,[0 25])
set(ax_in2,'xtick',-2:2)
set(ax_in2,'xticklabelrotation',0)
grid(ax_in2,'on')

% panel 3: m_3eqn,coarse
plot(ax(3),[0 10],[0 10],'k-','linewidth',lw)
clear h3
errorbar(ax(3),[mod_berg.m_mean(:,5); mod_glac.m_mean(:,5)],m(inc),m_ci_l(inc),m_ci_h(inc),'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
h3(4) = scatter(ax(3),mod_berg.m_mean(:,5),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
for i = 1:size(idx_glac,2)
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    h3(i) = scatter(ax(3),mod_glac.m_mean(idx_glac(find(any(idx_glac,2),1):end,i),5),m(idx_glac(:,i)),ms,dep_clrs(i,:),smbs{i},'filled','markeredgecolor','k');
end
% histogram
ax_in3 = axes(fig,'position',[ax(3).Position(1:2)+hist_loc 0 0]+[0 0 hist_size],'units','normalized','fontsize',fs-3,'xaxislocation','top');
hold(ax_in3,'on')
box(ax_in3,'on')
histogram(ax_in3,log10(r_crs),hist_lvls,'facecolor',0.*[1 1 1],'facealpha',1)
histogram(ax_in3,log10(r_crs(any(idx_glac(inc,:),2))),hist_lvls,'facecolor',0.7*[1 1 1],'facealpha',1)
% plot(ax_in3,log10(median(r,'omitnan'))*[1 1],[0 15],'k-','linewidth',1)
xlim(ax_in3,[-0.5 2.2])
ylim(ax_in3,[0 25])
set(ax_in3,'xtick',-2:2)
set(ax_in3,'xticklabelrotation',0)
grid(ax_in3,'on')

% axes limits
linkaxes(ax(2:3),'x')
xlim(ax(2),[0 4.25])
xlim(ax(1),[0 0.32])

linkaxes(ax,'y')
ylim(ax(1),[0 4.25])

% axes labels
ylabel(ax(1),'m_{obs} [m/day]','fontsize',fs)
xlabel(ax(1),'$\overline{|\vec{u}|}$ [m/s]','fontsize',fs,'interpreter','latex')
xlabel(ax(2),'m_{3eqn} [m/day]','fontsize',fs)
xlabel(ax(3),'m_{3eqn,coarse} [m/day]','fontsize',fs)
xlabel(ax_in2,'log_{10}(r_i)','fontsize',fs-3,'fontweight','normal')
ylabel(ax_in2,'count','fontsize',fs-3)
xlabel(ax_in3,'log_{10}(r_i)','fontsize',fs-3,'fontweight','normal')
ylabel(ax_in3,'count','fontsize',fs-3)

% legend
lgd = legend(ax(3),h3,cat(2,dep_lbls,{'berg'}),'orientation','horizontal','fontsize',fs-2);
lgd.Position([1 2]) = [ax(2).Position(1) sum(ax(2).Position([2 4]))+.01];

% panel lbls
for i = 1:length(ax)
    text(ax(i),.02,1,[char(96+i) ')'],'units','normalized','fontsize',fs-1,'horizontalalignment','left','verticalalignment','top')
end

% line labels
for i = 2:3
    text(ax(i),0.9,0.85,'1:1','units','normalized','fontsize',fs-2,'rotation',45,'horizontalalignment','center','verticalalignment','middle')
end
text(ax(2),0.65,0.9,'1.5:1','units','normalized','fontsize',fs-2,'rotation',atand(1.5),'horizontalalignment','center','verticalalignment','middle')

