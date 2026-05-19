% Script to plot simultaneous velocities measured by the meltstake ADCP and
% ADV.
%
% KJW
% 9 Mar 2025

clear

raw_dir = 'F:meltstake/data/raw';
proc_dir = 'F:meltstake/data/proc';
ms_tbl = loadMSInfo('segments');

% choose deployment and segment(s)
dep_num = 28;
seg_num = [1 Inf];

% find row number and time window
row_num = find((ms_tbl.Number == dep_num) & (ms_tbl.Window >= seg_num(1) & ms_tbl.Window <= seg_num(2)));
t1 = ms_tbl.Start(row_num(1));
t2 = ms_tbl.End(row_num(end));

% load
dep_name = ms_tbl.Folder{row_num(1)};
load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'))
load(fullfile(raw_dir,dep_name,'adv','adv.mat'))
load(fullfile(proc_dir,dep_name,'adv','svol_in_ocean.mat'))

% time indexing
idxt_adcp = adcp.burst.time>=t1 & adcp.burst.time<=t2;
idxt_adv = adv.time>=t1 & adv.time<=t2;

% transform ADV velocity
adv = msADVTransform(adv,adcp.attitude);

% fix ADCP coordinate system to be consistent with Weiss et al 2025 (+y
% outward from ice)
% adcp.burst.vel_ice(:,:,1) = -adcp.burst.vel_ice(:,:,1); % u
% adcp.burst.vel_ice(:,:,3) = -adcp.burst.vel_ice(:,:,3); % v

% ADV sample location (might need to fudge a little if ADCP is
% contaminated)
y_adv = 0.35; % m
[~,r_idx] = min(abs(adcp.burst.range-y_adv));

% extract ADCP velocity at ADV range
vel_adcp = adcp.burst.vel_ice(idxt_adcp,r_idx,[1 3 2]);
vel_adcp_xyz = adcp.burst.vel_xyz(idxt_adcp,r_idx,[1 2 5]); % [v_x v_y v_z]

% despike a little bit (bubbles...)
for i = 1:3
    vel_adcp(:,i) = medianFilter(vel_adcp(:,i),3);
    vel_adcp_xyz(:,i) = medianFilter(vel_adcp_xyz(:,i),3);
end

% ADV instrument vel corresponding to ADCP instrument coordinates
vel_adv_xyz = [-adv.vel_xyz(:,2) adv.vel_xyz(:,1) adv.vel_xyz(:,3)];
vel_adv_xyz(~idx_ocean,:) = nan;

%% plot
% plot to check for phase wrapping
figure(1); clf
plot(adv.time(idxt_adv),adv.vel(idxt_adv,:),'.-')
adcpQuickPlot(figure(2),adcp,'vel',extrema(adcp.burst.vel),[t1 t2],[0 1],1)

% plot ADCP velocity for reference
adcp.burst.vel_ice_uvw = adcp.burst.vel_ice(:,:,[1 3 2]);
adcpQuickPlot(figure(3),adcp,'vel_ice_uvw',0.2,[t1 t2],[0 0.65],1)

%% plot ADV and ADCP point velocity (ice coords)
figure(4); clf
pad = [.05 .05 .1 .1];
clear ax1
for i = 3:-1:1
    ax1(i) = axes(figure(4),'position',axgridpos(3,1,i,pad),'fontsize',10);
    hold on
    plot(adv.time(idxt_adv),adv.vel_ice(idxt_adv,i),'.-','color',colors(1))
    plot(adcp.burst.time(idxt_adcp),vel_adcp(:,i),'k-','linewidth',1)
    box on
    grid on
    text(0.02,0.93,sprintf('%c',char(116+i)),'units','normalized','fontsize',12)
    ylabel('vel [m/s]')
end
legend(ax1(1),{'ADV',sprintf('ADCP (y=%0.2f m)',round(y_adv,2))},'location','best','fontsize',10)
title(sprintf('deployment %d (%s)',dep_num,strrep(dep_name,'_','\_')),'fontsize',11)
set(ax1(1:2),'xticklabels',[])

linkaxes(ax1,'x')
xlim(ax1(1),[t1 t2])

%% plot ADV and ADCP point velocity (adcp inst coords)
figure(5); clf
pad = [.05 .05 .1 .1];
clear ax2
for i = 3:-1:1
    ax2(i) = axes(figure(5),'position',axgridpos(3,1,i,pad),'fontsize',10);
    hold on
    plot(adv.time(idxt_adv),vel_adv_xyz(idxt_adv,i),'.-','color',colors(1))
    plot(adcp.burst.time(idxt_adcp),vel_adcp_xyz(:,i),'k-','linewidth',1)
    box on
    grid on
    text(0.02,0.93,sprintf('v_%c',char(119+i)),'units','normalized','fontsize',12)
    ylabel('vel [m/s]')
end
legend(ax2(1),{'ADV',sprintf('ADCP (y=%0.2f m)',round(y_adv,2))},'location','best','fontsize',10)
title(sprintf('deployment %d (%s)',dep_num,strrep(dep_name,'_','\_')),'fontsize',11)
set(ax2(1:2),'xticklabels',[])

linkaxes(ax2,'x')
xlim(ax2(1),[t1 t2])
