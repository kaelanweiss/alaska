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
idx_berg = ms_tbl.Number <= 19;
for i = 1:3
    idx_glac(:,i) = ms_tbl.Number == (i+25);
end

% load model results
mod_glac = load('../../../data/melt_models_glacier.mat');
mod_berg = load('../../../data/melt_models.mat');

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

pad = [.035 .0 .06 .12];
shift = [0.35*pad(3) 0.6*pad(4)];

fig = figure(1); clf
setFigureSize(fig,[15 7])

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
% cbar.Label.String = 'T [K]';
clim(ax(1),[5 6.5])
text(ax(1),0.24,0.9,'(T-T_f) [K]','fontsize',fs-2,'horizontalalignment','center')

% panel 2: m_3eqn
plot(ax(2),[0 10],[0 10],'k-','linewidth',lw)
errorbar(ax(2),[mod_berg.m_mean(:,1); mod_glac.m_mean(:,1)],m,m_ci_l,m_ci_h,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
scatter(ax(2),mod_berg.m_mean(:,1),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
for i = 1:size(idx_glac,2)
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    scatter(ax(2),mod_glac.m_mean(idx_glac(find(any(idx_glac,2),1):end,i),1),m(idx_glac(:,i)),ms,clrs(i,:),smbs{i},'filled','markeredgecolor','k');
end

% panel 3: m_3eqn,coarse
plot(ax(3),[0 10],[0 10],'k-','linewidth',lw)
clear h3
errorbar(ax(3),[mod_berg.m_mean(:,5); mod_glac.m_mean(:,5)],m,m_ci_l,m_ci_h,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
h3(4) = scatter(ax(3),mod_berg.m_mean(:,5),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
for i = 1:size(idx_glac,2)
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    h3(i) = scatter(ax(3),mod_glac.m_mean(idx_glac(find(any(idx_glac,2),1):end,i),5),m(idx_glac(:,i)),ms,clrs(i,:),smbs{i},'filled','markeredgecolor','k');
end

% axes limits
linkaxes(ax(2:3),'x')
xlim(ax(2),[0 2.5])
xlim(ax(1),[0 0.32])

linkaxes(ax,'y')
ylim(ax(1),[0 4.25])

% axes labels
for i = 2:3
%     set(ax(i),'yticklabels',{})
end
ylabel(ax(1),'m_{obs} [m/day]','fontsize',fs)
xlabel(ax(1),'$\overline{\mathrm{U}_\mathrm{max}}$ [m/s]','fontsize',fs,'interpreter','latex')
xlabel(ax(2),'m_{3eqn} [m/day]','fontsize',fs)
xlabel(ax(3),'m_{3eqn,coarse} [m/day]','fontsize',fs)

% legend
lgd = legend(ax(3),h3,cat(2,dep_lbls,{'berg'}),'location','southeast','fontsize',fs-2);
lgd.Position = lgd.Position + [.035 -.015 0 0];
