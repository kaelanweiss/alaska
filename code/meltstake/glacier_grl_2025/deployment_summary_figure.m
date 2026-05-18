% Script to plot up Jonathan's idea for the full deployment time series.
%
% KJW
% 28 Jan 2025
clear

% raw dir
raw_dir = 'F:/meltstake/data/raw';

% choose deployment number
seg_tbl = loadMSInfo(28,'segments');
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

%% plot
fig = figure(1); clf
setFigureSize(fig,0.8*[21 24]);
clear ax

lw = 1;
fs = 9;

pad = [.05 .025 .15 .05];
shift = [-.08 0];
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
hold(ax(4),'on')
h4(4) = plot(ax(4),adcp.burst.time,hannFilter(umax,8*2),'-','color',0.2*[1 1 1],'linewidth',lw);
% for i = 1:3
%     h4(i) = plot(ax(4),adcp.burst.time,uvw_mean(:,i),'-','color',vel_clrs(i,:),'linewidth',lw);
% end

% T and S
hold(ax(5),'on')
for i = 1:3
    h5(i) = plot(ax(5),T(i).time,T(i).values,'-','color',T_clrs(4-i,:));
end
yyaxis(ax(5),'right')
set(ax(5),'ycolor','k')
h5(4) = plot(ax(5),hobo.time,hobo.S,'-','color',0.7*colors(6),'linewidth',1.5*lw);
yyaxis(ax(5),'left')

% melt
hold(ax(6),'on')
% 3eqn
h6(3) = plot(ax(6),t_3eqn,hannFilter(m_3eqn,2*2),'-','color',[1 0.6 0.6],'linewidth',lw);
% segments
for i = 1:nseg
    if i == 1
        h6(1) = plot([seg_tbl.Start(i) seg_tbl.End(i)],m(i)*[1 1],'-','linewidth',1*lw,'color',colors(1));
    else
        plot([seg_tbl.Start(i) seg_tbl.End(i)],m(i)*[1 1],'-','linewidth',1*lw,'color',colors(1))
    end
end
% mean 3eqn
h6(4) = plot(ax(6),t_melt+minutes(0.),m_3eqn_mean,'ko','markersize',6,'markerfacecolor','r','linewidth',lw);
% observed melt
errorbar(ax(6),t_melt,m,m_ci_l,m_ci_h,'.','linewidth',1.5*lw,'color','k')
h6(2) = plot(ax(6),t_melt,m,'ko','markersize',4,'linewidth',lw,'markerfacecolor','k');

% x tick
for i = 2:length(ax)-1
    set(ax(i),'xticklabel',{});
end

% ylabels
for i = 1:3
    ylabel(ax(i),'y [m]','fontsize',fs)
end
ylabel(ax(4),'[m/s]','fontsize',fs)
ylabel(ax(5),'T [\circC]','fontsize',fs)
yyaxis(ax(5),'right'); ylabel(ax(5),'S [psu]','fontsize',fs); yyaxis(ax(5),'left')
ylabel(ax(6),'melt [m/day]','fontsize',fs)

% panels
p_lbls = addPanelLabels(ax,fs+1);
for i = 1:length(p_lbls)
    % p_lbls(i).BackgroundColor = 'w';
    p_lbls(i).Position(1) = .01;
    p_lbls(i).Position(2) = 1.15;
    p_lbls(i).FontWeight = 'bold';
end
p_lbls(1).String = 'u'; p_lbls(2).String = 'v'; p_lbls(3).String = 'w';
% p_lbls(4).String = 'spatial mean vel';
p_lbls(4).String = 'sqrt(u^2+w^2)';
p_lbls(5).String = 'T and S';
p_lbls(6).String = 'melt rate';

% colobar
cbar = colorbar(ax(3),'position',cbarpos(ax(1:3),.02,.02));
cbar.Label.String = '[m/s]';
cbar.Label.FontSize = fs;

% legends
% lgd4 = legend(ax(4),h4,{'u','v','w','|U|_{max}'},'fontsize',fs);
% lgd4.Position(1) = 0.83;

lgd5 = legend(ax(5),h5,{'T(10cm)','T(25cm)','T(40cm)','S(30cm)'},'fontsize',fs);
lgd5.Position(1) = 0.83;

lgd6 = legend(ax(6),h6,{'segment','observed','3eqn','mean(3eqn)'},'fontsize',fs);
lgd6.Position(1) = 0.8;

% axes limits
linkaxes(ax,'x')
xlim(ax(1),[seg_tbl.Start(1) seg_tbl.End(end)]+minutes([-1 1]))

for i = 1:3
    ylim(ax(i),[0 y(1)])
    clim(ax(i),0.4*[-1 1])
end

% ylim(ax(4),0.55*[-1 1])

% ylim(ax(5),[3 6])
yyaxis(ax(5),'right')
% ylim(ax(5),[25 27])
yyaxis(ax(5),'left')
ylim(ax(6),[0 max(m+m_ci_h,[],'omitnan')+0.1])