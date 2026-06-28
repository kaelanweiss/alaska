% Script to plot up meltstake data for figure 2
%
% KJW
% 26 Jun 2026
clear

raw_dir = 'F:/meltstake/data/raw';
proc_dir = 'F:/meltstake/data/proc';

% choose deployment number
seg_tbl = loadMSInfo(28,'segments');
seg_tbl = seg_tbl(seg_tbl.include,:);
depname = seg_tbl.Folder{1};
nseg = size(seg_tbl,1);

% load colors
load glacier_clrs.mat

% load melt rate
m = msTable2Vector(seg_tbl.m)*0.24;
[m_ci_l,m_ci_h] = msTable2Vector(seg_tbl.m_ci);
m_ci_l = 0.24*m_ci_l; m_ci_h = 0.24*m_ci_h;
t_melt = mean([seg_tbl.Start seg_tbl.End],2);

% load T and S and adcp
load(fullfile(raw_dir,depname,'rbr','T.mat'))
load(fullfile(raw_dir,depname,'hobo','hobo.mat'))
hobo = hobo(2);
load(fullfile(raw_dir,depname,'adcp','adcp.mat'))

% load ctd
load(fullfile(proc_dir,'glacier_2024_ctd.mat'))
ctd = ctd(seg_tbl.Number(1)-25);

%% delete a little bit of bad S data
bad_sal_times = [datetime(2024,7,17,0,58,58.5) datetime(2024,7,17,1,1,7.5);...
                 datetime(2024,7,17,2,22,18.5) datetime(2024,7,17,2,22,36.5)];
idx_bad_sal = false(size(hobo.time));
for i = 1:size(bad_sal_times,1)
    t1 = bad_sal_times(i,1);
    t2 = bad_sal_times(i,2);
    idx_bad_sal = idx_bad_sal | (hobo.time>=t1 & hobo.time<=t2);
end
% hobo.S(idx_bad_sal) = interp1(hobo.time(~idx_bad_sal),hobo.S(~idx_bad_sal),hobo.time(idx_bad_sal),'linear');
hobo.S(idx_bad_sal) = nan;

%% calculate spatial mean velocities
adcp.burst.range = cast(adcp.burst.range,'double');
ridx = adcp.burst.range <= min(seg_tbl.rmax-0.03);
y = mean(seg_tbl.rmax)/cosd(25)-adcp.burst.range(ridx);
adcp = msADCPTransform(adcp,adcp.burst.processing.cor_min,adcp.burst.processing.amp_min);


uvw = adcp.burst.vel_ice(:,ridx,[1 3 2]);
uvw_mean = squeeze(mean(uvw,2,'omitnan'));
for i = 1:size(uvw_mean,2)
    uvw_mean(:,i) = hannFilter(uvw_mean(:,i),8*2);
end

uvw_filt = uvw;
for i = 1:size(uvw,1)
    for k = 1:size(uvw,3)
        uvw_filt(i,:,k) = hannFilter(uvw(i,:,k),3);
    end
end
umax = max(vecnorm(uvw_filt(:,:,[1 3]),2,3),[],2);

%% plot 1
fig = figure(1); clf
setFigureSize(fig,[17.5 14]);
clear ax

lw = 1;
fs = 11;

pad = [.05 .025 .1 .05];
shift = [-.022 0.03];
nf = 6;
for i = 1:nf
    ax(i) = axes(fig,'position',axgridpos(nf,1,i,pad,shift),'fontsize',fs);
    box(ax(i),'on')
end

% vel pcolors
uw = uvw(:,:,[1 3]);
for i = 1:2
    pcolor(ax(i),adcp.burst.time,y,uw(:,:,i)')
    shading(ax(i),'flat')
    cmocean('bal',ax(i))
    clim(ax(i),0.4*[-1 1])
end

% vel magnitude
plot(ax(3),adcp.burst.time,hannFilter(umax,8*2),'-','color',0.2*[1 1 1],'linewidth',lw)
set(ax(3),'fontsize',fs)

% T
hold(ax(4),'on')
plot(ax(4),extrema(T(1).time),ctd.T_mean(ctd.zidx)*[1 1],'color',0.6*[1 1 1],'LineWidth',3*lw)
% plotConfidenceRegion(ax(4),extrema(T(1).time),T_outer(1)*[1 1],T_outer(2)*[1 1],'k',0.3);
for i = 1:3
    plot(ax(4),T(i).time,T(i).values,'-','color',T_clrs(4-i,:),'linewidth',lw);
end
for i = 1:3
    h4(i) = plot(ax(4),extrema(T(i).time),[0 0],'-','color',T_clrs(4-i,:),'linewidth',2*lw);
end

% S
hold(ax(5),'on')
plot(ax(5),extrema(T(1).time),ctd.S_mean(ctd.zidx)*[1 1],'color',0.6*[1 1 1],'LineWidth',3*lw)
% plotConfidenceRegion(ax(5),extrema(hobo.time),S_outer(1)*[1 1],S_outer(2)*[1 1],0.0*colors(6),0.3);
plot(ax(5),hobo.time,hobo.S,'-','color',0.7*colors(6),'linewidth',1.5*lw);

% melt
hold(ax(6),'on')
for i = 1:nseg
    if i == 1
        h6(1) = plot([seg_tbl.Start(i) seg_tbl.End(i)],m(i)*[1 1],'-','linewidth',1*lw,'color',0.5*[1 1 1]);
    else
        plot([seg_tbl.Start(i) seg_tbl.End(i)],m(i)*[1 1],'-','linewidth',1*lw,'color',0.5*[1 1 1])
    end
end
errorbar(ax(6),t_melt,m,m_ci_l,m_ci_h,'.','linewidth',1.5*lw,'color','k')
h6(2) = plot(ax(6),t_melt,m,'ko','markersize',4,'linewidth',lw,'markerfacecolor','k');

% % umax
% yyaxis(ax(6),'right')
% set(ax(6),'ycolor','k')
% plot(ax(6),adcp.burst.time,hannFilter(umax,8*5),'r')
% ylim(ax(6),[0 0.45])
% yyaxis(ax(6),'left')

% x tick
for i = 1:length(ax)-1
    set(ax(i),'xticklabel',{});
end

% ylabels
for i = 1:2
    ylabel(ax(i),'y [m]','fontsize',fs)
    set(ax(i),'fontsize',fs,'ytick',0:0.2:0.6)
end
ylabel(ax(3),'speed [m/s]','fontsize',fs)
ylabel(ax(4),'T [\circC]','fontsize',fs)
ylabel(ax(5),'S [psu]','fontsize',fs)
ylabel(ax(6),'m [m/day]','fontsize',fs)

% panels
p_lbls = addPanelLabels(ax,fs);
for i = 1:length(p_lbls)
    % p_lbls(i).BackgroundColor = 'w';
    p_lbls(i).Position(1) = .005;
    p_lbls(i).Position(2) = 1.02;
    p_lbls(i).FontWeight = 'bold';
end
p_lbls(1).String = 'a) u';
p_lbls(2).String = 'b) w';
p_lbls(3).String = 'c) sqrt(u^2+w^2)';
p_lbls(4).String = 'd) temperature';
p_lbls(5).String = 'e) salinity';
p_lbls(6).String = 'f) melt rate';

% colobar
cbar = colorbar(ax(2),'position',cbarpos(ax(1:2),.01,.02));
cbar.Label.String = 'velocity [m/s]';
cbar.Label.FontSize = fs;
% cbar.Label.FontWeight = 'bold';
cbar.FontSize = fs-1;

% legends
lgd4 = legend(ax(4),h4,{'10cm','25cm','40cm'},'fontsize',fs-1,'orientation','horizontal','location','southeast');
lgd4.Position(1:2) = [0.5 0.385];

% axes limits
linkaxes(ax,'x')
xlim(ax(1),[seg_tbl.Start(1) seg_tbl.End(end)]+0.1*minutes([-1 1]))

for i = 1:2
    ylim(ax(i),[0 y(1)])
    clim(ax(i),0.4*[-1 1])
end

ylim(ax(4),[5 7.5])
ylim(ax(5),[25.5 28])
ylim(ax(6),[0 4.5])
