% Script to plot simultaneous velocities measured by the meltstake ADCP and
% ADV.
%
% KJW
% 9 Mar 2025

clear

raw_dir = 'F:meltstake/data/raw';

tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

% choose deployment and segment
dep_num = 27;
seg_num = 1;

% find row number and time window
row_num = (ms_tbl.Number == dep_num) & (ms_tbl.Window == seg_num);
t1 = ms_tbl.Start(row_num);
t2 = ms_tbl.End(row_num);

% load
dep_path = fullfile(raw_dir,ms_tbl.Folder{row_num});
load(fullfile(dep_path,'adcp','adcp.mat'))
load(fullfile(dep_path,'adv','adv.mat'))

% time indexing
idxt_adcp = adcp.burst.time>=t1 & adcp.burst.time<=t2;
idxt_adv = adv.time>=t1 & adv.time<=t2;

% transform ADV velocity (until my processing catches up)
attitude = adcp.attitude;
% pitch = attitude.roll-270;
% roll = -attitude.pitch;
% attitude.pitch = pitch;
% attitude.roll = roll;
adv = msADVTransform(adv,attitude);

% fix ADCP coordinate system to be consistent with Weiss et al 2025 (+y
% outward from ice)
adcp.burst.vel_ice(:,:,1) = -adcp.burst.vel_ice(:,:,1); % u
adcp.burst.vel_ice(:,:,3) = -adcp.burst.vel_ice(:,:,3); % v

% ADV sample location (might need to fudge a little if ADCP is
% contaminated)
y_adv = 0.40; % m
[~,r_idx] = min(abs(adcp.burst.range-y_adv));

% extract ADCP velocity at ADV range
vel_adcp = adcp.burst.vel_ice(idxt_adcp,r_idx,[1 3 2]);
% despike a little bit (bubbles...)
for i = 1:3
    vel_adcp(:,i) = medianFilter(vel_adcp(:,i),3);
end
%% plot
% plot to check for phase wrapping
figure(1); clf
plot(adv.time(idxt_adv),adv.vel(idxt_adv,:),'.-')
adcpQuickPlot(figure(2),adcp,'vel',extrema(adcp.burst.vel),[t1 t2],[0 1],1)

% plot ADCP velocity for reference
adcp.burst.vel_ice_uvw = adcp.burst.vel_ice(:,:,[1 3 2]);
adcpQuickPlot(figure(3),adcp,'vel_ice_uvw',0.2,[t1 t2],[0 0.65],1)

%% plot ADV and ADCP point velocity 
figure(4); clf
pad = [.05 .05 .1 .1];
clear ax
for i = 3:-1:1
    ax(i) = axes(figure(4),'position',axgridpos(3,1,i,pad),'fontsize',10);
    hold on
    plot(adv.time(idxt_adv),adv.vel_ice(idxt_adv,i),'.-','color',colors(1))
    plot(adcp.burst.time(idxt_adcp),vel_adcp(:,i),'k-','linewidth',1)
    box on
    grid on
    text(0.02,0.93,sprintf('%c',char(116+i)),'units','normalized','fontsize',12)
    ylabel('vel [m/s]')
end
legend(ax(1),{'ADV',sprintf('ADCP (y=%0.2f m)',round(y_adv,2))},'location','best','fontsize',10)
title(sprintf('deployment %d (%s) segment %d',dep_num,strrep(ms_tbl.Folder{row_num},'_','\_'),seg_num),'fontsize',11)
set(ax(1:2),'xticklabels',[])

linkaxes(ax,'x')
xlim(ax(1),[t1 t2])
