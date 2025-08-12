% Script to create figure 3 for the melt rate paper. This is a 20 tile 
% layout comparing the T, U, u, w, and v of 4 different sections of 
% meltstake data.
%
% KJW
% 4 Sep 2024

set(0,'defaulttextinterpreter','tex');

clear

addpath('..')

% load data
proc_path = 'F:meltstake\data\proc';
ms_tbl = loadMSInfo(1:19,'manualwindows');
load ../../../data/melt_models.mat

% define sections by deployment index (table row)
sections = [8 3 29 9];
section_lbls = {'crossflow','plume','wave','eddy'};

% zoom into time limits
t0 = [datetime(2023,7,8,22,50,0);...
      datetime(2023,5,29,21,18,0);...
      datetime(2023,9,23,21,38,0);...
      datetime(2023,7,9,1,22,0)];

t_width = [minutes(5);...
           minutes(5);...
           minutes(2);...
           minutes(5)];

% melt rate
m = msTable2Vector(ms_tbl.m);
m_ci = ms_tbl.m_ci;
m = m(sections)*.24;
m_ci = m_ci(sections)*.24;
m_3eqn = m_mean(sections,1);

m_err = m./m_3eqn;
m_err_ci = m_ci./m_3eqn;

% preallocate data space
ns = length(sections);
T_all = cell(ns,1);
adcp_all = cell(ns,1);
rmax_all = 0;

% hanning width
hann_dt = 2; % s

for i = 1:ns    
    % load section data
    dep_num = ms_tbl.Number(sections(i));
    wind_num = ms_tbl.Window(sections(i));
    fprintf('%d.%d\n',dep_num,wind_num)
    
    load(fullfile(proc_path,ms_tbl.Folder{sections(i)},sprintf('T%d.mat',wind_num)))
    load(fullfile(proc_path,ms_tbl.Folder{sections(i)},sprintf('adcp%d.mat',wind_num)))

    adcp = msADCPTransform(adcp,adcp.burst.processing.cor_min,adcp.burst.processing.amp_min);

    % flip order of thermistors if need be
    if any(sections(i)==[15])
        T.T = fliplr(T.T);
    end
    
    % convert ADCP range to distance from ice (estimate)
    rmax = ms_tbl.rmax(sections(i));
    adcp.burst.range = rmax/cosd(25)-adcp.burst.range;% + 0.5*adcp.burst.cellsize;
    rmax_all = max([rmax_all adcp.burst.range(1)]);

    % calculate velocity scale
    vel = adcp.burst.vel_ice(:,:,1:2); % u,w
    % smooth across 3 points in the profile
    for j = 1:size(vel,1)
        for p = 1:size(vel,3)
            vel(j,:,p) = hannFilter(squeeze(vel(j,:,p)),3);
        end
    end
    % calculate magnitude
    vel_mag = vecnorm(vel,2,3);
    % profile max and mean
    vel_max = max(vel_mag,[],2);
    adcp.burst.vel_scale = vel_max;

    % calculate velocity profiles
    adcp.burst.vel_prof = squeeze(mean(cat(3,vel,vel_mag),1,'omitnan'));
    adcp.burst.vel_prof_std = squeeze(std(cat(3,vel,vel_mag),1,'omitnan'));

    % smoothing
    k_T = round(hann_dt/seconds(diff(T.time(1:2))));
    k_u = hann_dt*adcp.burst.samplerate;
    % T
    for j = 1:size(T.T,2)
        T.T(:,j) = hannFilter(T.T(:,j),k_T);
    end
    % U
    adcp.burst.vel_scale = hannFilter(adcp.burst.vel_scale,k_u);
    for j = 1:size(adcp.burst.vel,2)
        for p = 1:2
            adcp.burst.vel_ice(:,j,p) = hannFilter(adcp.burst.vel_ice(:,j,p),k_u);
        end
    end

    T_all{i} = T;
    adcp_all{i} = adcp;

end

%% plot
% plot params
% padding
pad = [.01 .025 .07 .11];
shift = [-.01 0];

lw = 1.2;
fs = 11;

% T colors
T_clr = [0.2 1 0.5]/1.2;
T_clrs = T_clr'*[1 0.5 0];
T_lgd_lbls = {'near','mid','far'};

% labels
panel_lbls = {'T_w [\circC]','U_{max} [m/s]','u [m/s]','w [m/s]','v [m/s]'};
% panel_lbls = {'$T_w$','max($|\vec{u}|$)','$u$','$w$','$v$'};
y_lbls = {'[\circC]','[m/s]',{'y [m]'},'y [m]',{'y [m]'}};
% y_lbls = {'$^\circ$C','m/s',{''},'y, distance from ice (m)',{''}};

% manual lims
T_lim = [3 9];
vel_lim = [0 0.24];
pc_lims = [.15 .15 .15 .15]-.049;

% create and size figure
figsize = [16 16];
fig = figure(4); clf
setFigureSize(fig,figsize);
clear ax

% loop through each section
ncols = 5;
for i = 1:ns
    % plot
    % T
    ax(ncols*(i-1)+1) = axes(fig,'position',axgridpos(ns,5,ncols*(i-1)+1,pad,shift));
    hold on
    for j = 1:3
        plot(T_all{i}.T(:,j),minutes(T_all{i}.time-t0(i)),'linewidth',lw,'color',T_clrs(:,j))
    end
    xlim(T_lim)
    box on
    ax(ncols*(i-1)+1).YTick = 0:5;

    % U_max
    adcp_time = minutes(adcp_all{i}.burst.time-t0(i));
    ax(ncols*(i-1)+2) = axes(fig,'position',axgridpos(ns,5,ncols*(i-1)+2,pad,shift));
    plot(hannFilter(adcp_all{i}.burst.vel_scale,8),adcp_time,'k','linewidth',lw)
    xlim(vel_lim)

    % u
    ax(ncols*(i-1)+3) = axes(fig,'position',axgridpos(ns,5,ncols*(i-1)+3,pad,shift));
    pcolor(adcp_all{i}.burst.range,adcp_time,adcp_all{i}.burst.vel_ice(:,:,1))
    shading flat
    cmocean('bal')
    clim(pc_lims(i)*[-1 1])
    xlim([0 rmax_all])

    % w
    ax(ncols*(i-1)+4) = axes(fig,'position',axgridpos(ns,5,ncols*(i-1)+4,pad,shift));
    pcolor(adcp_all{i}.burst.range,adcp_time,adcp_all{i}.burst.vel_ice(:,:,2))
    shading flat
    cmocean('bal')
    clim(pc_lims(i)*[-1 1])
    xlim([0 rmax_all])

    % v
    ax(ncols*(i-1)+5) = axes(fig,'position',axgridpos(ns,5,ncols*(i-1)+5,pad,shift));
    pcolor(adcp_all{i}.burst.range,adcp_time,adcp_all{i}.burst.vel_ice(:,:,3))
    shading flat
    cmocean('bal')
    clim(pc_lims(i)*[-1 1])
    xlim([0 rmax_all])

    % link axes
    linkaxes(ax((1:ncols)+ncols*(i-1)),'y')
    ylim(ax(ncols*i),[0 minutes(t_width(i))])

end

set(ax,'FontSize',fs-2)

% clear internal ticks
% x
set(ax(1:(ns-1)*ncols),'xticklabels',{})

% y
set(ax(setxor(1:(ns*ncols),1:ncols:ns*ncols)),'yticklabels',{})

% colorbars
cbar = colorbar(ax(ncols));
cbar.Position = cbarpos(ax(ncols),.005,.015);
cbar.FontSize = fs-2;
text(ax(5),1.16,1.1,'m/s','fontsize',fs-2,'units','normalized')

% panel and axis labels
patch_dx = .35;
patch_dy = .2;
for i = 1:ncols
    for j = 1:ns
        txt = text(ax(i+ncols*(j-1)),1-.03,.97,['(' char(96 + ncols*(j-1) + i) ')'],'units','normalized','fontsize',fs,...
        'verticalalignment','top','horizontalalignment','right');
        if i == 3 && j == 1
            txt.BackgroundColor = 'w';
        end
        if i == 1
            ylabel(ax(i+ncols*(j-1)),'time [min]','fontsize',fs)
        end
    end
    title(ax(i),panel_lbls{i},'fontsize',fs,'FontWeight','normal')
    xlabel(ax(i+(ns-1)*ncols),y_lbls{i},'fontsize',fs)
end

% T legend
axpos = ax(ncols+1).Position;
lgd = legend(ax(ncols+1),T_lgd_lbls,'orientation','vertical','location','southeast',...
    'fontsize',fs-3,'edgecolor','k');
ax(ncols+1).Position = axpos;
% lgd.Position = lgd.Position - [.05 -.005 0 0];

