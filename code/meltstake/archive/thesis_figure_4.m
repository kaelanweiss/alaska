% Script to create new figure 4 for the thesis/melt rate paper. This is a
% 20 tile layout comparing the T, U, u, w, and v of 4 different sections of
% meltstake data.
%
% KJW
% 11 Apr 2024
set(0,'defaulttextinterpreter','latex');

clear

% load data
proc_path = 'F:meltstake\data\proc';
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');
load ../../data/melt_models.mat

% define sections by deployment index (table row)
sections = [8 3 29 9];
section_lbls = {'crossflow','plume','wave','eddy'};

% [9 15;...
%          8 29;...
%          9 22;...
%          9 21;...
%          9 3;...
%          9 18];%...
%          25 31;...
%          26 31];

% zoom into time limits
t0 = [datetime(2023,7,8,22,50,0);...
      datetime(2023,5,29,21,18,0);...
      datetime(2023,9,23,21,38,0);...
      datetime(2023,7,9,1,22,0)];
%       datetime(2023,7,11,1,8,0);...
%       datetime(2023,7,10,19,30,0);...
%       datetime(2023,7,8,22,50,0) datetime(2023,9,23,21,38,0);...
%       datetime(2023,7,9,1,40,0) datetime(2023,7,11,1,4,0);...
%       datetime(2023,7,9,1,40,0) datetime(2023,7,11,1,8,0);...
%       datetime(2023,7,9,1,22,0) datetime(2023,5,29,21,18,0);...
%       NaT NaT];
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
pxi = .015;
pxo = .09;
pyi = .015;
pyo = .1;

lw = 1.2;
fs = 11;

% T colors
T_clr = [0.2 1 0.5]/1.2;
T_clrs = T_clr'*[1 0.5 0];
T_lgd_lbls = {'near','mid','far'};

% labels
panel_lbls = {'$T$','max($|\vec{u}|$)','$u$','$w$','$v$'};
y_lbls = {'$^\circ$C','m/s',{''},'y, distance from ice (m)',{''}};

% manual lims
T_lim = [3 9];
vel_lim = [0 0.24];
pc_lims = [.15 .15 .15 .15]-.049;

% create and size figure
figsize = [18 16];
fig = figure(4); clf
set(fig,'units','centimeters'); 
fpos = get(fig,'position');
fpos(3:4) = figsize;
set(fig,'position',fpos);
set(fig,'paperunits','centimeters')
set(fig,'papersize',figsize)
clear ax

% loop through each section
nrows = 5;
for i = 1:ns
    % plot
    % T
    ax(nrows*(i-1)+1) = axes(fig,'position',axgridpos(5,ns,nrows*(i-1)+1,pxi,pyi,pxo,pyo,'flip'));
    hold on
    for j = 1:3
        plot(minutes(T_all{i}.time-t0(i)),T_all{i}.T(:,j),'linewidth',lw,'color',T_clrs(:,j))
    end
    title({sprintf('%s (%s)',char(64+i),section_lbls{i}),...
        strrep(sprintf('$m_{obs}$: %.1f$\\pm$%.1f$\\frac{m}{day}$',m(i),m_ci(i)),'m$0.','m$.'), strrep(sprintf('$r_{3eqn}$: %.1f$\\pm$%.1f',m_err(i),m_err_ci(i)),'m$0.','m$.')},'fontsize',fs)
    ylim(T_lim)
    box on

    % vel scale
    adcp_time = minutes(adcp_all{i}.burst.time-t0(i));
    ax(nrows*(i-1)+2) = axes(fig,'position',axgridpos(5,ns,nrows*(i-1)+2,pxi,pyi,pxo,pyo,'flip'));
    plot(adcp_time,hannFilter(adcp_all{i}.burst.vel_scale,8),'k','linewidth',lw)
    ylim(vel_lim)

    % u
    ax(nrows*(i-1)+3) = axes(fig,'position',axgridpos(5,ns,nrows*(i-1)+3,pxi,pyi,pxo,pyo,'flip'));
    pcolor(adcp_time,adcp_all{i}.burst.range,-adcp_all{i}.burst.vel_ice(:,:,1)')
    shading flat
    cmocean('bal')
    clim(pc_lims(i)*[-1 1])
    ylim([0 rmax_all])

    % w
    ax(nrows*(i-1)+4) = axes(fig,'position',axgridpos(5,ns,nrows*(i-1)+4,pxi,pyi,pxo,pyo,'flip'));
    pcolor(adcp_time,adcp_all{i}.burst.range,adcp_all{i}.burst.vel_ice(:,:,2)')
    shading flat
    cmocean('bal')
    clim(pc_lims(i)*[-1 1])
    ylim([0 rmax_all])

    % w
    ax(nrows*(i-1)+5) = axes(fig,'position',axgridpos(5,ns,nrows*(i-1)+5,pxi,pyi,pxo,pyo,'flip'));
    pcolor(adcp_time,adcp_all{i}.burst.range,-adcp_all{i}.burst.vel_ice(:,:,3)')
    shading flat
    cmocean('bal')
    clim(pc_lims(i)*[-1 1])
    ylim([0 rmax_all])

    % link axes
    linkaxes(ax((1:nrows)+nrows*(i-1)),'x')
    xlim(ax(nrows*i),[0 minutes(t_width(i))])


end

% clear internal ticks
% x
set(ax(setxor(1:(nrows*ns),[nrows:nrows:(nrows*ns)])),'xticklabels',{})
% y
set(ax(setxor(1:(nrows*ns),1:nrows)),'yticklabels',{})

% colorbars
cbu = colorbar(ax(nrows*ns-2));
cbw = colorbar(ax(nrows*ns));
cbu.Position = cbarpos(ax(nrows*ns-2),.005,.015);
cbw.Position = cbarpos(ax(nrows*ns),.005,.015);
cbw.Position(4) = sum(cbu.Position([2 4])) - cbw.Position(2);
delete(cbu);
cbw.Label.String = 'm/s';
cbw.Label.Interpreter = 'latex';
cbw.Label.FontSize = fs;

% panel and axis labels
patch_dx = .35;
patch_dy = .2;
for i = 1:nrows
    for j = 1:ns
        txt = text(ax(i+nrows*(j-1)),.025,.97,['(' char(96 + nrows*(j-1) + i) ') ' panel_lbls{i}],'units','normalized','fontsize',fs,...
        'verticalalignment','top','horizontalalignment','left');
        if i >= 3
            txt.BackgroundColor = 'w';
        end
        if i == nrows
            xlabel(ax(i+nrows*(j-1)),'time (min)','fontsize',fs)
        end
    end
    ylabel(ax(i),y_lbls{i},'fontsize',fs)
end

% ax_xy_lbl = axes(fig,'position',axgridpos(2,1,2,pxi,pyi,pxo,pyo));
% xlabel(ax_xy_lbl,'time (min)','fontsize',fs)
% ylabel(ax_xy_lbl,'distance from ice (m)','fontsize',fs)
% set(ax_xy_lbl,'off')

% T legend
axpos = ax(nrows*(ns-1)+1).Position;
lgd = legend(ax(nrows*(ns-1)+1),T_lgd_lbls,'orientation','vertical','location','southeastoutside',...
    'fontsize',fs-3,'edgecolor','k');
ax(nrows*(ns-1)+1).Position = axpos;
lgd.Position = lgd.Position - [.05 -.005 0 0];

