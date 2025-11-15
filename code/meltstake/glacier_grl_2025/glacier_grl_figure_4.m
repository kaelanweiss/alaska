% Create a plot to compare momentum, temperature, salinity, and melt rate
% between Polly, meltstakes, and moorings(?)
%
% KJW
% 10 Nov 2025
clear

addpath('..')

% load data
load platform_comparison_data\polly_data.mat
seg_tbl = loadMSInfo(26:28,'segments');

dep_nums = unique(seg_tbl.Number);
ndeps = length(dep_nums);
dep_idx = false(size(seg_tbl,1),ndeps);
for i = 1:ndeps
    dep_idx(:,i) = seg_tbl.Number==dep_nums(i);
end

u_ms = msTable2Vector(seg_tbl.u_max);
T_ms = msTable2Vector(seg_tbl.T);
S_ms = msTable2Vector(seg_tbl.S);
m_ms = msTable2Vector(seg_tbl.m)*0.24;
d_ms = seg_tbl.depth;

%% plot
lw = 1;
fs = 11;
msize = 5;
mstyle = 'x';
alpha = 0.6;

xlbls = {'ocean speed [m/s]','temperature [\circC]','salinity [psu]','melt rate [m/day]'};

pad = [.02 .0 .06 .12];
shift = [0.35*pad(3) 0.6*pad(4)];

fig = figure(1); clf
setFigureSize(fig,[15 7]);

clear ax h
for i = 1:4
    ax(i) = axes(fig,'position',axgridpos(1,4,i,pad,shift),'fontsize',fs-2);
    box on
    hold on
    axis ij
end

% a) momentum
plotConfidenceRegion(ax(1),adcp_depth,mom_mean-mom_std,mom_mean+mom_std,colors(1),alpha,'x');
h(1) = plot(ax(1),mom_mean,adcp_depth,'color',colors(1),'linewidth',lw);
h(2) = plot(ax(1),u_ms,round(d_ms),'k','markersize',msize,'marker',mstyle,'linestyle','none','linewidth',lw);

% b) temperature
plotConfidenceRegion(ax(2),ctd_depth,T_mean(:,1),T_mean(:,3),colors(2),alpha,'x');
plot(ax(2),T_mean(:,2),ctd_depth,'color',colors(2),'linewidth',lw)
plot(ax(2),T_ms,round(d_ms),'k','markersize',msize,'marker',mstyle,'linestyle','none','linewidth',lw)

% c) salinity
plotConfidenceRegion(ax(3),ctd_depth,S_mean(:,1),S_mean(:,3),colors(3),alpha,'x');
plot(ax(3),S_mean(:,2),ctd_depth,'color',colors(3),'linewidth',lw)
plot(ax(3),S_ms,round(d_ms),'k','markersize',msize,'marker',mstyle,'linestyle','none','linewidth',lw)

% d) melt rate
plot(ax(4),m_3eqn*86400,adcp_depth,'k-','linewidth',lw)
plot(ax(4),2*m_3eqn*86400,adcp_depth,'k--','linewidth',lw)
plot(ax(4),m_ms,round(d_ms),'k','markersize',msize,'marker',mstyle,'linestyle','none','linewidth',lw)


% limits
linkaxes(ax,'y')
ylim(ax(1),[0 120])
xlim(ax(1),[0 0.48])
xlim(ax(2),[4 8])
xlim(ax(3),[24.8 29.2])
xlim(ax(4),[0 5.5])

% labels
for i = 2:4
    set(ax(i),'yticklabels','')
end
ylabel(ax(1),'depth [m]','fontsize',fs)
for i = 1:4
    xlabel(ax(i),xlbls{i},'fontsize',fs)
end

% legend
lgd = legend(ax(2),h,{'RHIB','MS'},'fontsize',fs-2,'location','southeast');
lgd.Position([1 2]) = [0.2 0.22];
