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
dep_num = 27;
load(fullfile(proc_dir,'glacier_2024_ctd.mat'));
ctd_all = ctd;
ctd = ctd_all(dep_num-25);
idx_depth = ctd.depth_mean >= 1;

% load segment T and S
seg_tbl = loadMSInfo(dep_num,'segments');
T = msTable2Vector(seg_tbl.T);
S = msTable2Vector(seg_tbl.S);
z_ms = mean(seg_tbl.depth);

%% plot
fig = figure(1); clf
setFigureSize(fig,[6 5]);
clear ax

lw = 1.8;
fs = 10;
alpha = 0.5;
ms = 3;
S_clr = 0.7*[35 142 35]/255; % 0.7*colors(6)

pad = [.03 .025 .11 .15];
shift = [0.095 0.082];
nf = 2;
for i = 1:nf
    ax(i) = axes(fig,'position',axgridpos(1,nf,i,pad,shift),'fontsize',fs);
    box(ax(i),'on')
    hold(ax(i),'on')
    % plot(ax(i),[0 50],ctd.depth(ctd.zidx)*[1 1],'-','color',0.6*[1 1 1],'linewidth',2)
end

% temperature
c_num = 2;
plotConfidenceRegion(ax(1),ctd.depth_mean(idx_depth),ctd.T_min(idx_depth),ctd.T_max(idx_depth),T_clrs(c_num,:),alpha,'x');
% plot(ax(1),ctd.T,ctd.depth,'color',T_clrs(c_num,:),'linewidth',lw/2)
plot(ax(1),ctd.T_mean(idx_depth),ctd.depth_mean(idx_depth),'color',T_clrs(c_num,:),'linewidth',lw);
plot(ax(1),T,z_ms*ones(size(T)),'kdiamond','linewidth',0.5*lw,'markersize',ms);

% salinity
plotConfidenceRegion(ax(2),ctd.depth_mean(idx_depth),ctd.S_min(idx_depth),ctd.S_max(idx_depth),S_clr,alpha,'x');
plot(ax(2),ctd.S_mean(idx_depth),ctd.depth_mean(idx_depth),'color',S_clr,'linewidth',lw)
plot(ax(2),S,z_ms*ones(size(S)),'kdiamond','linewidth',0.5*lw,'markersize',ms)

% flip
for i = 1:nf
    axis(ax(i),'ij')
    set(ax(i),'fontsize',fs)
end

% axis limits
switch dep_num
    case 26
        xlim(ax(1),[3.5 7.5])
        xlim(ax(2),[23.8 29])
    case 27
        xlim(ax(1),[4 8])
        xlim(ax(2),[23.3 31.2])
    case 28
        xlim(ax(1),[4.5 7.8])
        xlim(ax(2),[22.5 32])
end
linkaxes(ax,'y')
ylim(ax(1),[0 150])

% ticks
set(ax(1),'xtick',0:10)
set(ax(2),'yticklabel',{})
set(ax(2),'xtick',24:2:30)

% axis labels
xlabel(ax(1),'T [\circC]','fontsize',fs)
xlabel(ax(2),'S [psu]','fontsize',fs)
ylabel(ax(1),'depth [m]','fontsize',fs)

% panel labels
p_lbls = addPanelLabels(ax,fs);
p_lbls(1).String = 'h) fjord T';
p_lbls(2).String = 'i) fjord S';
for i = 1:2
    p_lbls(i).FontWeight = 'bold';
    p_lbls(i).Position(1:2) = [.005 1.1];
end

% text
% ms_txt = text(ax(1),min(xlim(ax(1)))+0.1,ctd.depth(ctd.zidx),{'MS','depth'},'fontsize',fs,'verticalalignment','top');