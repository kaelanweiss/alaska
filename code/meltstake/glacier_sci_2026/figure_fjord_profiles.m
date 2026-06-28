% Script to plot open-fjord profiles during meltstake deployments
%
% KJW
% 26 Jun 2026

clear

addpath('..')

proc_dir = 'F:/meltstake/data/proc';

% load colors
load ../glacier_grl_2025/glacier_clrs.mat

% load ctd data
load(fullfile(proc_dir,'glacier_2024_ctd.mat'));
ctd_all = ctd;
ctd = ctd_all(3);

%% plot
fig = figure(1); clf
setFigureSize(fig,[7.5 8]);
clear ax

lw = 2;
fs = 11;
alpha = 0.7;

pad = [.05 .025 .1 .12];
shift = [0.07 0.07];
nf = 2;
for i = 1:nf
    ax(i) = axes(fig,'position',axgridpos(1,nf,i,pad,shift),'fontsize',fs);
    box(ax(i),'on')
    hold(ax(i),'on')
    plot(ax(i),[0 50],ctd.depth(ctd.zidx)*[1 1],'k--')
end

% temperature
c_num = 2;
plotConfidenceRegion(ax(1),ctd.depth_mean,ctd.T_min,ctd.T_max,T_clrs(c_num,:),alpha,'x');
% plot(ax(1),ctd.T,ctd.depth,'color',T_clrs(c_num,:),'linewidth',lw/2)
plot(ax(1),ctd.T_mean,ctd.depth_mean,'color',T_clrs(c_num,:),'linewidth',lw)

% salinity
plotConfidenceRegion(ax(2),ctd.depth_mean,ctd.S_min,ctd.S_max,0.7*colors(6),alpha,'x');
plot(ax(2),ctd.S_mean,ctd.depth_mean,'color',0.7*colors(6),'linewidth',lw)

% flip
for i = 1:nf
    axis(ax(i),'ij')
    set(ax(i),'fontsize',fs)
end

% axis limits
xlim(ax(1),[4.5 7.8])
xlim(ax(2),[22.5 32])
linkaxes(ax,'y')
ylim(ax(1),[0 150])

% ticks
set(ax(2),'yticklabel',{})
set(ax(2),'xtick',24:2:30)

% axis labels
xlabel(ax(1),'T [\circC]','fontsize',fs)
xlabel(ax(2),'S [psu]','fontsize',fs)
ylabel(ax(1),'depth [m]','fontsize',fs)

% panel labels
p_lbls = addPanelLabels(ax,fs);
p_lbls(1).String = 'g) fjord T';
p_lbls(2).String = 'h) fjord S';
for i = 1:2
    p_lbls(i).FontWeight = 'bold';
end

% text
ms_txt = text(ax(1),min(xlim(ax(1)))+0.1,ctd.depth(ctd.zidx),{'MS','depth'},'fontsize',fs,'verticalalignment','top');