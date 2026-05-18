% Script to plot up meltstake data 
%
% KJW
% 28 Jan 2025
clear

% raw dir
raw_dir = 'F:/meltstake/data/raw';

% choose deployment number
seg_tbl = loadMSInfo(28,'segments');
seg_tbl = seg_tbl(seg_tbl.include,:);
depname = seg_tbl.Folder{1};
nseg = size(seg_tbl,1);

% load colors
load ../glacier_grl_2025/glacier_clrs.mat

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

%% fix a little bit of S data
bad_sal_times = [datetime(2024,7,17,0,58,58.5) datetime(2024,7,17,1,1,7.5);...
                 datetime(2024,7,17,2,22,18.5) datetime(2024,7,17,2,22,36.5)];
idx_bad_sal = false(size(hobo.time));
for i = 1:size(bad_sal_times,1)
    t1 = bad_sal_times(i,1);
    t2 = bad_sal_times(i,2);
    idx_bad_sal = idx_bad_sal | (hobo.time>=t1 & hobo.time<=t2);
end
hobo.S(idx_bad_sal) = interp1(hobo.time(~idx_bad_sal),hobo.S(~idx_bad_sal),hobo.time(idx_bad_sal),'linear');


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

% calculate 3eqn melt
t_3eqn = T(1).time;
T_3eqn = T(1).values;
u_3eqn = interp1(adcp.burst.time,umax,t_3eqn);
S_3eqn = interp1(hobo.time,hobo.S,t_3eqn,'nearest',median(hobo.S,'omitnan'));
m_3eqn = solve3Eqn(u_3eqn,T_3eqn,S_3eqn,0)*86400;

m_3eqn_mean = nan(nseg,1);
m_3eqn_std = nan(nseg,1);
for i = 1:nseg
    idxi = t_3eqn>=seg_tbl.Start(i) & t_3eqn<=seg_tbl.End(i);
    m_3eqn_mean(i) = mean(m_3eqn(idxi),'omitnan');
    m_3eqn_std(i) = std(m_3eqn(idxi),'omitnan');
end

%% plot 1
fig = figure(1); clf
setFigureSize(fig,[22 18]);
clear ax

lw = 1;
fs = 11;

% outer T and S
T_outer = [6.8 7.2];
S_outer = [27.4 27.7];

pad = [.05 .025 .08 .05];
shift = [-.022 0.02];
nf = 6;
for i = 1:nf
    ax(i) = axes(fig,'position',axgridpos(nf,1,i,pad,shift),'fontsize',fs);
    box(ax(i),'on')
end

% vel pcolors
for i = 1:3
    pcolor(ax(i),adcp.burst.time,y,uvw(:,:,i)')
    shading(ax(i),'flat')
    cmocean('bal',ax(i))
    clim(ax(i),0.4*[-1 1])
end

% y-mean vel
% hold(ax(4),'on')
% h4(4) = plot(ax(4),adcp.burst.time,hannFilter(umax,8*2),'-','color',0.2*[1 1 1],'linewidth',lw);
% for i = 1:3
%     h4(i) = plot(ax(4),adcp.burst.time,uvw_mean(:,i),'-','color',vel_clrs(i,:),'linewidth',lw);
% end

% T
hold(ax(4),'on')
plotConfidenceRegion(ax(4),extrema(T(1).time),T_outer(1)*[1 1],T_outer(2)*[1 1],'k',0.3);
for i = 1:3
    plot(ax(4),T(i).time,T(i).values,'-','color',T_clrs(4-i,:),'linewidth',lw);
end
for i = 1:3
    h4(i) = plot(ax(4),extrema(T(i).time),[0 0],'-','color',T_clrs(4-i,:),'linewidth',2*lw);
end

% S
hold(ax(5),'on')
plotConfidenceRegion(ax(5),extrema(hobo.time),S_outer(1)*[1 1],S_outer(2)*[1 1],0.0*colors(6),0.3);
plot(ax(5),hobo.time,hobo.S,'-','color',0.7*colors(6),'linewidth',1.5*lw);

% melt
hold(ax(6),'on')
% 3eqn
% h6(3) = plot(ax(6),t_3eqn,hannFilter(m_3eqn,2*2),'-','color',[1 0.6 0.6],'linewidth',lw);
% segments
for i = 1:nseg
    if i == 1
        h6(1) = plot([seg_tbl.Start(i) seg_tbl.End(i)],m(i)*[1 1],'-','linewidth',1*lw,'color',0.5*[1 1 1]);
    else
        plot([seg_tbl.Start(i) seg_tbl.End(i)],m(i)*[1 1],'-','linewidth',1*lw,'color',0.5*[1 1 1])
    end
end
% mean 3eqn
% h6(4) = plot(ax(6),t_melt+minutes(0.),m_3eqn_mean,'ko','markersize',6,'markerfacecolor','r','linewidth',lw);
% observed melt
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
for i = 1:3
    ylabel(ax(i),'y [m]','fontsize',fs,'fontweight','bold')
    set(ax(i),'fontsize',fs,'ytick',0:0.2:0.6)
end
ylabel(ax(4),'T [\circC]','fontsize',fs,'fontweight','bold')
ylabel(ax(5),'S [psu]','fontsize',fs,'fontweight','bold')
ylabel(ax(6),'melt [m/day]','fontsize',fs,'fontweight','bold')

% panels
p_lbls = addPanelLabels(ax,fs+1);
for i = 1:length(p_lbls)
    % p_lbls(i).BackgroundColor = 'w';
    p_lbls(i).Position(1) = .01;
    p_lbls(i).Position(2) = 1.15;
    p_lbls(i).FontWeight = 'bold';
end
p_lbls(1).String = 'u (horz. along-ice)';
p_lbls(2).String = 'v (horz. ice-normal)';
p_lbls(3).String = 'w (vert.)';
p_lbls(4).String = 'temperature';
p_lbls(5).String = 'salinity';
p_lbls(6).String = 'melt rate';

% colobar
cbar = colorbar(ax(3),'position',cbarpos(ax(1:3),.02,.02));
cbar.Label.String = '[m/s]';
cbar.Label.FontSize = fs;
cbar.Label.FontWeight = 'bold';
cbar.FontSize = fs-1;

% legends

lgd4 = legend(ax(4),h4,{'10cm','25cm','40cm'},'fontsize',fs,'orientation','horizontal','location','southeast');
lgd5.Position(1) = 0.83;

% lgd6 = legend(ax(6),h6,{'segment','observed','3eqn','mean(3eqn)'},'fontsize',fs);
% lgd6.Position(1) = 0.8;

% axes limits
linkaxes(ax,'x')
xlim(ax(1),[seg_tbl.Start(1) seg_tbl.End(end)]+0.1*minutes([-1 1]))

for i = 1:3
    ylim(ax(i),[0 y(1)])
    clim(ax(i),0.4*[-1 1])
end

ylim(ax(4),[5 7.5])
ylim(ax(5),[25 28])
ylim(ax(6),[0 max(m+m_ci_h,[],'omitnan')+0.1])

%% plot 2
fig2 = figure(2); clf
setFigureSize(fig2,[26 20]);
clear ax

lw = 1.5;
fs = 12;

% outer T and S
T_outer = [6.8 7.2];
S_outer = [27.4 27.7];

pad = [.05 .03 .08 .05];
shift = [-.022 0.02];
nf = 6;
for i = 1:nf
    ax(i) = axes(fig2,'position',axgridpos(nf,1,i,pad,shift),'fontsize',fs);
    box(ax(i),'on')
end

% melt
hold(ax(1),'on')

m_3eqn_clr = 0.9*vel_clrs(3,:);
% 3eqn time series
% plot(ax(1),t_3eqn,hannFilter(m_3eqn,3),'-','linewidth',lw,'color',m_3eqn_clr+0.2*([1 1 1]-m_3eqn_clr));

% 3eqn
% h6(3) = plot(ax(6),t_3eqn,hannFilter(m_3eqn,2*2),'-','color',[1 0.6 0.6],'linewidth',lw);
% segments
for i = 1:nseg
    if i == 1
        h6(1) = plot(ax(1),[seg_tbl.Start(i) seg_tbl.End(i)]+0.7*minutes([1 -1]),0.25*[1 1],'|-','linewidth',1*lw,'color',0.3*[1 1 1]);
    else
        plot(ax(1),[seg_tbl.Start(i) seg_tbl.End(i)]+0.7*minutes([1 -1]),0.25*[1 1],'|-','linewidth',1*lw,'color',0.3*[1 1 1])
    end
end
% mean 3eqn
% h6(4) = plot(ax(6),t_melt+minutes(0.),m_3eqn_mean,'ko','markersize',6,'markerfacecolor','r','linewidth',lw);
% observed melt
errorbar(ax(1),t_melt,m,m_ci_l,m_ci_h,'.','linewidth',1.2*lw,'color','k','capsize',10)
h6(2) = plot(ax(1),t_melt,m,'k^','markersize',8,'linewidth',0.5*lw,'markerfacecolor',dep_clrs(3,:));

% T
hold(ax(2),'on')
plotConfidenceRegion(ax(2),extrema(T(1).time),T_outer(1)*[1 1],T_outer(2)*[1 1],T_clrs(2,:),0.6);
for i = 1:3
    plot(ax(2),T(i).time,T(i).values,'-','color',T_clrs(4-i,:),'linewidth',lw);
end
for i = 1:3
    h4(i) = plot(ax(2),extrema(T(i).time),[0 0],'-','color',T_clrs(4-i,:),'linewidth',2*lw);
end
% plot(ax(2),T(1).time,T(2).values,'-','color',T_clrs(4-3,:),'linewidth',lw)

% S
hold(ax(3),'on')
plotConfidenceRegion(ax(3),extrema(hobo.time),S_outer(1)*[1 1],S_outer(2)*[1 1],0.7*vel_clrs(3,:),0.6);
plot(ax(3),hobo.time,hobo.S,'-','color',0.7*vel_clrs(3,:),'linewidth',1.2*lw);

% vel mag
plot(ax(4),adcp.burst.time,hannFilter(umax,8*2),'-','color',vel_clrs(1,:),'linewidth',lw)
set(ax(4),'fontsize',fs)

% vel pcolors
pcolor(ax(5),adcp.burst.time,y,uvw(:,:,1)')
pcolor(ax(6),adcp.burst.time,y,uvw(:,:,3)')
for i = 5:6
    shading(ax(i),'flat')
    cmocean('bal',ax(i))
    clim(ax(i),0.4*[-1 1])
end

% y-mean vel
% hold(ax(4),'on')
% h4(4) = plot(ax(4),adcp.burst.time,hannFilter(umax,8*2),'-','color',0.2*[1 1 1],'linewidth',lw);
% for i = 1:3
%     h4(i) = plot(ax(4),adcp.burst.time,uvw_mean(:,i),'-','color',vel_clrs(i,:),'linewidth',lw);
% end

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
for i = 5:6
    set(ax(i),'fontsize',fs,'ytick',0:0.2:0.6)
    ylabel(ax(i),'y [m]','fontsize',fs,'fontweight','bold')
end
ylabel(ax(2),'T [\circC]','fontsize',fs,'fontweight','bold')
ylabel(ax(3),'S [psu]','fontsize',fs,'fontweight','bold')
ylabel(ax(4),'speed [m/s]','fontsize',fs,'fontweight','bold')
ylabel(ax(1),'melt [m/day]','fontsize',fs,'fontweight','bold')

% panels
p_lbls = addPanelLabels(ax,fs+3);
for i = 1:length(p_lbls)
    % p_lbls(i).BackgroundColor = 'w';
    p_lbls(i).Position(1) = .01;
    p_lbls(i).Position(2) = 1.2;
    p_lbls(i).FontWeight = 'bold';
end
p_lbls(4).Position(2) = 1.22;
p_lbls(4).String = 'sqrt(u^2+w^2)';
p_lbls(5).String = 'u (horizontal along-ice)';
% p_lbls(5).String = 'v (horz. ice-normal)';
p_lbls(6).String = 'w (vertical)';
p_lbls(2).String = 'temperature';
p_lbls(3).String = 'salinity';
p_lbls(1).String = 'melt rate';

% colobar
cbar = colorbar(ax(6),'position',cbarpos(ax(5:6),.005,.015));
cbar.Label.String = '[m/s]';
cbar.Label.FontSize = fs;
cbar.Label.FontWeight = 'bold';
cbar.FontSize = fs-1;
cbar.FontWeight = 'bold';

% legends
% lgd4 = legend(ax(2),h4,{'10cm','25cm','40cm'},'fontsize',fs,'orientation','horizontal','location','southeast');
lgd4 = legend(ax(2),h4,{'10cm','25cm','40cm'},'fontsize',fs,'fontweight','bold');
lgd4.Position(1) = sum([ax(2).Position([1 3]) .005]);
lgd4.Position(2) = get(h4(1),'parent').Position(2);

% lgd6 = legend(ax(6),h6,{'segment','observed','3eqn','mean(3eqn)'},'fontsize',fs);
% lgd6.Position(1) = 0.8;

% axes limits
linkaxes(ax,'x')
xlim(ax(1),[seg_tbl.Start(1) seg_tbl.End(end)]+0.1*minutes([-1 1]))

for i = 4:6
    ylim(ax(i),[0 y(1)])
    clim(ax(i),0.4*[-1 1])
end

ylim(ax(2),[5 7.5])
ylim(ax(3),[25.5 28.2])
ylim(ax(4),[0 0.48])
ylim(ax(1),[0 max(m+m_ci_h,[],'omitnan')+0.1])

%% melt rate timeseries
fig3 = figure(3); clf
setFigureSize(fig3,[20 4.5]);

lw = 1.5;
fs = 12;

ax3 = axes(fig3,'position',axgridpos(1,1,1,[0 0 .15 .15],[-0.08 0.09]),'fontsize',fs);

hold(ax3,'on')
box(ax3,'on')

m_3eqn_clr = 0.9*vel_clrs(3,:);

% 3eqn time series
h3(2) = plot(t_3eqn,hannFilter(m_3eqn,3),'-','linewidth',lw,'color',m_3eqn_clr+0.5*([1 1 1]-m_3eqn_clr));
xlim(ax3,[seg_tbl.Start(1) seg_tbl.End(end)]+0.1*minutes([-1 1]))
ylabel(ax3,{'melt rate [m/day]'},'fontsize',fs,'fontweight','bold')

% segments
for i = 1:nseg
    if i == 1
        h3(1) = plot(ax3,[seg_tbl.Start(i) seg_tbl.End(i)]+minutes([2 -2]),0.2*[1 1],'|-','linewidth',1*lw,'color',0.3*[1 1 1]);
    else
        plot(ax3,[seg_tbl.Start(i) seg_tbl.End(i)]+0.7*minutes([1 -1]),0.2*[1 1],'|-','linewidth',1*lw,'color',0.3*[1 1 1])
    end
end

% observed melt
h3(5) = errorbar(ax3,t_melt,m,m_ci_l,m_ci_h,'.','linewidth',1.5*lw,'color','k');
h3(4) = plot(ax3,t_melt,m,'k^','markersize',7,'linewidth',0.8*lw,'markerfacecolor',dep_clrs(3,:));

% 3 eqn means
h3(3) = plot(ax3,t_melt+minutes(0.),m_3eqn_mean,'khexagram','markersize',8,'markerfacecolor',m_3eqn_clr,'linewidth',0.8*lw);


lgd3 = legend(ax3,h3(1:4),{'segment','model','mean(model)','observed'},'fontsize',fs,'fontweight','bold');
lgd3.Position(1) = 0.79;

ylim(ax3,[0 4.1])