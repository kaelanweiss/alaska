% Melt rate plot for the glacier Science paper. Modified from the earlier
% GRL version.
%
% KJW
% 19 Jun 2026

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
mod_glac = load('../../../data/melt_models_glacier_hann_1sec.mat');
mod_berg = load('../../../data/melt_models.mat');
mod_glac.m_mean = mod_glac.m_mean(any(idx_glac(ms_tbl.Start>datetime(2024,7,1),:),2),:);
mod_berg.m_mean = mod_berg.m_mean(idx_berg(ms_tbl.Start<datetime(2024,7,1)),:);

% cm/hr to m/dy
cmh2md = 0.24;

% extract fields
[m,~] = msTable2Vector(ms_tbl.m);
[T,T_ci] = msTable2Vector(ms_tbl.T);
Tb = ms_tbl.Tb;
[S,S_ci] = msTable2Vector(ms_tbl.S);
[u,u_ci] = msTable2Vector(ms_tbl.u_max);
[m_ci_l,m_ci_h] = msTable2Vector(ms_tbl.m_ci);

% convert melt rate to m/day
m = m*cmh2md;
m_ci_l = m_ci_l*cmh2md;
m_ci_h = m_ci_h*cmh2md;

% calculate transfer coefficient enhancement
coeff = mean(mod_glac.m_mean(:,10)./mod_glac.m_mean(:,1));

% uncertainty on tuning
tune_unc = .12;

% melt ratio
r = m(inc)./[coeff*mod_berg.m_mean(:,1); mod_glac.m_mean(:,10)];
r_ci = m_ci_l(inc)./[coeff*mod_berg.m_mean(:,1); mod_glac.m_mean(:,10)];
r_crs = m(inc)./[coeff*mod_berg.m_mean(:,5); mod_glac.m_mean(:,11)];
r_crs_ci = m_ci_l(inc)./[coeff*mod_berg.m_mean(:,5); mod_glac.m_mean(:,11)];

%% plot new
lw = 0.5;
fs = 10;
marker_scale = 30;
berg_marker = 'diamond';

hist_lvls = -1.5:0.15:2.3;
hist_loc = [.48 .00];
hist_size = [0.25 0.13];

pad = [.035 .08 .275 .05];
shift = [-0.6*pad(3) .03];

fig3 = figure(3); clf
setFigureSize(fig3,[7.5 14])

clear ax
for i = 1:3
    ax(i) = axes(fig3,'position',axgridpos(3,1,i,pad,shift),'fontsize',fs-2);
    box on
    hold on
end

dep_lbls = cell(1,3);
for i = 1:3
    dep_lbls{i} = sprintf('%dm',round(mean(ms_tbl.depth(idx_glac(:,i)))));
end

% panel 1: velocity and temperature
clear h1
errorbar(ax(1),u,m,m_ci_l,m_ci_h,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
scatter(ax(1),u(idx_berg),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
for i = 1:size(idx_glac,2)
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    h1(i) = scatter(ax(1),u(idx_glac(:,i)),m(idx_glac(:,i)),ms,T(idx_glac(:,i))-Tb(idx_glac(:,i)),smbs{i},'filled','markeredgecolor','k');
end
cmocean('sol',ax(1))
cbar = colorbar(ax(1),'position',cbarpos(ax(1),.01,.03),'fontsize',ax(1).FontSize);
clim(ax(1),[5.5 7])
text(ax(1),1.3,0.5,{'T-T_b','[K]'},'fontsize',fs-1,'horizontalalignment','center','verticalalignment','middle','units','normalized','fontweight','normal')

% panel 2: m_3eqn
% plotConfidenceRegion(ax(2),[0 10],[0 10]*(1-tune_unc),[0 10]*(1+tune_unc),'k',0.15)
plot(ax(2),[0 10],[0 10],'k-','linewidth',lw)
plot(ax(2),[0 10],(1/coeff)*[0 10],'k--','linewidth',lw)
errorbar(ax(2),[coeff*mod_berg.m_mean(:,1); mod_glac.m_mean(:,10)],m(inc),m_ci_l(inc),m_ci_h(inc),'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
scatter(ax(2),coeff*mod_berg.m_mean(:,1),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
for i = 1:size(idx_glac,2)
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    scatter(ax(2),mod_glac.m_mean(idx_glac(find(any(idx_glac,2),1):end,i),10),m(idx_glac(:,i)),ms,dep_clrs(i,:),smbs{i},'filled','markeredgecolor','k');
end
% histogram
ax_in2 = axes(fig3,'position',[ax(2).Position(1:2)+hist_loc 0 0]+[0 0 hist_size],'units','normalized','fontsize',fs-2,'yaxislocation','right');
hold(ax_in2,'on')
box(ax_in2,'on')
histogram(ax_in2,log10(r),hist_lvls,'facecolor',0.7*[1 1 1],'facealpha',1)
histogram(ax_in2,log10(r(any(idx_glac(inc,:),2))),hist_lvls,'facecolor',0.*[1 1 1],'facealpha',1)
% plot(ax_in2,log10(median(r,'omitnan'))*[1 1],[0 15],'k-','linewidth',1)
xlim(ax_in2,[-0.5 2.2])
ylim(ax_in2,[0 25])
set(ax_in2,'xtick',-2:2)
set(ax_in2,'xticklabelrotation',0)
grid(ax_in2,'on')

% panel 3: m_3eqn,coarse
% plotConfidenceRegion(ax(3),[0 10],[0 10]*(1-tune_unc),[0 10]*(1+tune_unc),'k',0.15)
plot(ax(3),[0 10],[0 10],'k-','linewidth',lw)
clear h3
errorbar(ax(3),[coeff*mod_berg.m_mean(:,5); mod_glac.m_mean(:,11)],m(inc),m_ci_l(inc),m_ci_h(inc),'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
h3(4) = scatter(ax(3),coeff*mod_berg.m_mean(:,5),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
for i = 1:size(idx_glac,2)
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    h3(i) = scatter(ax(3),mod_glac.m_mean(idx_glac(find(any(idx_glac,2),1):end,i),11),m(idx_glac(:,i)),ms,dep_clrs(i,:),smbs{i},'filled','markeredgecolor','k');
end
% histogram
ax_in3 = axes(fig3,'position',[ax(3).Position(1:2)+hist_loc 0 0]+[0 0 hist_size],'units','normalized','fontsize',fs-2,'yaxislocation','right');
hold(ax_in3,'on')
box(ax_in3,'on')
histogram(ax_in3,log10(r_crs),hist_lvls,'facecolor',0.7*[1 1 1],'facealpha',1)
histogram(ax_in3,log10(r_crs(any(idx_glac(inc,:),2))),hist_lvls,'facecolor',0.*[1 1 1],'facealpha',1)
% plot(ax_in3,log10(median(r,'omitnan'))*[1 1],[0 15],'k-','linewidth',1)
xlim(ax_in3,[-0.5 2.2])
ylim(ax_in3,[0 25])
set(ax_in3,'xtick',-2:2)
set(ax_in3,'xticklabelrotation',0)
grid(ax_in3,'on')

% axes limits
linkaxes(ax(2:3),'x')
xlim(ax(2),[0 4.5])
xlim(ax(1),[0 0.32])

linkaxes(ax,'y')
ylim(ax(1),[0 4.6])

% axes labels
for i = 1:3
    ylabel(ax(i),'m_{obs} [m day^{-1}]','fontsize',fs)
end
xlabel(ax(1),'U_{inst} [m s^{-1}]','fontsize',fs)
xlabel(ax(2),'m_{3eqn} [m day^{-1}]','fontsize',fs)
xlabel(ax(3),'m_{3eqn} [m day^{-1}]','fontsize',fs)
xlabel(ax_in2,'log_{10}(m_{obs}/m_{3eqn})','fontsize',fs-1,'fontweight','normal')
ylabel(ax_in2,'segments','fontsize',fs-1)
xlabel(ax_in3,'log_{10}(m_{obs}/m_{3eqn})','fontsize',fs-1,'fontweight','normal')
ylabel(ax_in3,'segments','fontsize',fs-1)

% legend
lgd = legend(ax(3),h3,cat(2,dep_lbls,{'iceberg'}),'orientation','vertical','fontsize',fs-1);
lgd.Position([1 2]) = [sum(ax(2).Position([1 3]))+.01 sum(ax(2).Position([2 4]))-lgd.Position(4)/2];

% panel lbls
pnl_lbls = addPanelLabels([ax ax_in2 ax_in3],fs-1);
pnl_lbls(1).String = 'a) U_{inst} and T-T_b';
pnl_lbls(2).String = 'b) inst forcing';
pnl_lbls(3).String = 'd) mean forcing';
pnl_lbls(4).String = 'c)';
pnl_lbls(5).String = 'e)';
for i = 1:length(pnl_lbls)
    pnl_lbls(i).FontWeight = 'bold';
    pnl_lbls(i).Position(2) = 0.975;
end

% ticks
for i = 2:3
    set(ax(i),'xtick',0:4)
end

% line labels
for i = 2:3
    text(ax(i),0.9,0.82,'1:1','units','normalized','fontsize',fs-2,'rotation',45,'horizontalalignment','center','verticalalignment','middle')
end
text(ax(2),0.75,0.35,{'Jenkins','et al. (2010)'},'units','normalized','fontsize',fs-2,'rotation',atand(1/coeff),'horizontalalignment','center','verticalalignment','middle')

% rmspe labels
title(ax_in2,'RMSPE: 0.27','fontsize',fs-1,'fontweight','normal')
title(ax_in3,'RMSPE: 0.71','fontsize',fs-1,'fontweight','normal')


% %% plot
% lw = 0.5;
% fs = 10;
% marker_scale = 30;
% berg_marker = 'diamond';
% 
% hist_lvls = -1.5:0.15:2.3;
% hist_loc = [.48 .00];
% hist_size = [0.25 0.13];
% 
% pad = [.035 .08 .275 .05];
% shift = [-0.6*pad(3) .03];
% 
% fig2 = figure(2); clf
% setFigureSize(fig2,[7.5 14])
% 
% clear ax
% for i = 1:3
%     ax(i) = axes(fig2,'position',axgridpos(3,1,i,pad,shift),'fontsize',fs-2);
%     box on
%     hold on
% end
% 
% dep_lbls = cell(1,3);
% for i = 1:3
%     dep_lbls{i} = sprintf('%dm',round(mean(ms_tbl.depth(idx_glac(:,i)))));
% end
% 
% % panel 1: velocity and temperature
% clear h1
% errorbar(ax(1),u,m,m_ci_l,m_ci_h,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
% scatter(ax(1),u(idx_berg),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
% for i = 1:size(idx_glac,2)
%     ms = marker_scale;
%     if strcmp(smbs{i},'square')
%         ms = ms*1.3;
%     end
%     h1(i) = scatter(ax(1),u(idx_glac(:,i)),m(idx_glac(:,i)),ms,T(idx_glac(:,i))-Tb(idx_glac(:,i)),smbs{i},'filled','markeredgecolor','k');
% end
% cmocean('sol',ax(1))
% cbar = colorbar(ax(1),'position',cbarpos(ax(1),.01,.03),'fontsize',ax(1).FontSize);
% clim(ax(1),[5.5 7])
% text(ax(1),1.3,0.5,{'T-T_b','[K]'},'fontsize',fs-1,'horizontalalignment','center','verticalalignment','middle','units','normalized','fontweight','normal')
% 
% % panel 2: m_3eqn
% plot(ax(2),[0 10],[0 10],'k-','linewidth',lw)
% plot(ax(2),[0 10],(1/1.5)*[0 10],'k--','linewidth',lw)
% errorbar(ax(2),coeff*[mod_berg.m_mean(:,1); mod_glac.m_mean(:,1)],m(inc),m_ci_l(inc),m_ci_h(inc),'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
% scatter(ax(2),coeff*mod_berg.m_mean(:,1),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
% for i = 1:size(idx_glac,2)
%     ms = marker_scale;
%     if strcmp(smbs{i},'square')
%         ms = ms*1.3;
%     end
%     scatter(ax(2),coeff*mod_glac.m_mean(idx_glac(find(any(idx_glac,2),1):end,i),1),m(idx_glac(:,i)),ms,dep_clrs(i,:),smbs{i},'filled','markeredgecolor','k');
% end
% % histogram
% ax_in2 = axes(fig2,'position',[ax(2).Position(1:2)+hist_loc 0 0]+[0 0 hist_size],'units','normalized','fontsize',fs-2,'yaxislocation','right');
% hold(ax_in2,'on')
% box(ax_in2,'on')
% histogram(ax_in2,log10(r),hist_lvls,'facecolor',0.7*[1 1 1],'facealpha',1)
% histogram(ax_in2,log10(r(any(idx_glac(inc,:),2))),hist_lvls,'facecolor',0.*[1 1 1],'facealpha',1)
% % plot(ax_in2,log10(median(r,'omitnan'))*[1 1],[0 15],'k-','linewidth',1)
% xlim(ax_in2,[-0.5 2.2])
% ylim(ax_in2,[0 25])
% set(ax_in2,'xtick',-2:2)
% set(ax_in2,'xticklabelrotation',0)
% grid(ax_in2,'on')
% 
% % panel 3: m_3eqn,coarse
% plot(ax(3),[0 10],[0 10],'k-','linewidth',lw)
% clear h3
% errorbar(ax(3),coeff*[mod_berg.m_mean(:,5); mod_glac.m_mean(:,5)],m(inc),m_ci_l(inc),m_ci_h(inc),'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
% h3(4) = scatter(ax(3),coeff*mod_berg.m_mean(:,5),m(idx_berg),marker_scale/2,0.4*[1 1 1],berg_marker,'filled','markeredgecolor','k');
% for i = 1:size(idx_glac,2)
%     ms = marker_scale;
%     if strcmp(smbs{i},'square')
%         ms = ms*1.3;
%     end
%     h3(i) = scatter(ax(3),coeff*mod_glac.m_mean(idx_glac(find(any(idx_glac,2),1):end,i),5),m(idx_glac(:,i)),ms,dep_clrs(i,:),smbs{i},'filled','markeredgecolor','k');
% end
% % histogram
% ax_in3 = axes(fig2,'position',[ax(3).Position(1:2)+hist_loc 0 0]+[0 0 hist_size],'units','normalized','fontsize',fs-2,'yaxislocation','right');
% hold(ax_in3,'on')
% box(ax_in3,'on')
% histogram(ax_in3,log10(r_crs),hist_lvls,'facecolor',0.7*[1 1 1],'facealpha',1)
% histogram(ax_in3,log10(r_crs(any(idx_glac(inc,:),2))),hist_lvls,'facecolor',0.*[1 1 1],'facealpha',1)
% % plot(ax_in3,log10(median(r,'omitnan'))*[1 1],[0 15],'k-','linewidth',1)
% xlim(ax_in3,[-0.5 2.2])
% ylim(ax_in3,[0 25])
% set(ax_in3,'xtick',-2:2)
% set(ax_in3,'xticklabelrotation',0)
% grid(ax_in3,'on')
% 
% % axes limits
% linkaxes(ax(2:3),'x')
% xlim(ax(2),[0 4.5])
% xlim(ax(1),[0 0.32])
% 
% linkaxes(ax,'y')
% ylim(ax(1),[0 4.6])
% 
% % axes labels
% for i = 1:3
%     ylabel(ax(i),'m_{obs} [m day^{-1}]','fontsize',fs)
% end
% xlabel(ax(1),'U_{inst} [m s^{-1}]','fontsize',fs)
% xlabel(ax(2),'m_{3eqn} [m day^{-1}]','fontsize',fs)
% xlabel(ax(3),'m_{3eqn} [m day^{-1}]','fontsize',fs)
% xlabel(ax_in2,'log_{10}(m_{obs}/m_{3eqn})','fontsize',fs-1,'fontweight','normal')
% ylabel(ax_in2,'segments','fontsize',fs-1)
% xlabel(ax_in3,'log_{10}(m_{obs}/m_{3eqn})','fontsize',fs-1,'fontweight','normal')
% ylabel(ax_in3,'segments','fontsize',fs-1)
% 
% % legend
% lgd = legend(ax(3),h3,cat(2,dep_lbls,{'iceberg'}),'orientation','vertical','fontsize',fs-1);
% lgd.Position([1 2]) = [sum(ax(2).Position([1 3]))+.01 sum(ax(2).Position([2 4]))-lgd.Position(4)/2];
% 
% % panel lbls
% pnl_lbls = addPanelLabels([ax ax_in2 ax_in3],fs-1);
% pnl_lbls(1).String = 'a) U_{inst} and T-T_b';
% pnl_lbls(2).String = 'b) full forcing';
% pnl_lbls(3).String = 'd) mean forcing';
% pnl_lbls(4).String = 'c)';
% pnl_lbls(5).String = 'e)';
% for i = 1:length(pnl_lbls)
%     pnl_lbls(i).FontWeight = 'bold';
%     pnl_lbls(i).Position(2) = 0.975;
% end
% 
% % ticks
% for i = 2:3
%     set(ax(i),'xtick',0:4)
% end
% 
% % line labels
% for i = 2:3
%     text(ax(i),0.9,0.82,'1:1','units','normalized','fontsize',fs-2,'rotation',45,'horizontalalignment','center','verticalalignment','middle')
% end
% % text(ax(2),0.63,0.9,'1.6:1','units','normalized','fontsize',fs-2,'rotation',atand(1.5),'horizontalalignment','center','verticalalignment','middle')
% % text(ax(2),0.83,0.61,{'Jenkins','et al. (2010)'},'units','normalized','fontsize',fs-3,'rotation',atand(1/1.6),'horizontalalignment','center','verticalalignment','middle')
% text(ax(2),0.75,0.35,{'Jenkins','et al. (2010)'},'units','normalized','fontsize',fs-2,'rotation',atand(1/1.6),'horizontalalignment','center','verticalalignment','middle')

