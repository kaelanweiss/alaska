% Script to plot velocity from ADV within a deployment normalized by the
% ADCP "outer" velocity.
%
% notes: probably want to process adv into segments at some point and then
% update, also need to try to "unwind" velocity field
%
% KJW
% 22 May 2026
clear

proc_dir = 'F:/meltstake/data/proc';
raw_dir = 'F:/meltstake/data/raw';

% deployment number
dep_num = 28;
ms_tbl = loadMSInfo(dep_num,'segments');
nsegs = size(ms_tbl,1);
dep_name = ms_tbl.Folder{1};

% set adcp distance for outer velocity calculation
idxr = adcp.burst.range <= 0.2;

% load data
load(fullfile(proc_dir,dep_name,'adv','pck_edges.mat'))
load(fullfile(proc_dir,dep_name,'adv','svol_in_ocean.mat'))
load(fullfile(raw_dir,dep_name,'adv','adv.mat'),'adv')
load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'))
int_vel = load(fullfile(proc_dir,dep_name,'velocity_lowpass.mat'));

% transform velocities (may be redundant)
adcp = msADCPTransform(adcp,adcp.burst.processing.cor_min,adcp.burst.processing.amp_min);
adv = msADVTransform(adv,adcp.attitude);
% add adv magnitude and bin by burst
adv.vel_ice = [adv.vel_ice vecnorm(adv.vel_ice,2,2)];
adv = burstAverageADV(adv);
% length of burst
b_dt = adv.samples_per_burst/adv.sample_rate;

%% calculate velocity things
% set up arrays
y_adv = [];
u_norm = [];
w_norm = [];
u_outer = [];

y_max = 0;
for i = 1:nsegs
    fprintf('segment %d\n',i)
    % load adcp
    load(fullfile(proc_dir,dep_name,sprintf('adcp%d.mat',i)))
    % create indexing vectors
    t1 = ms_tbl.Start(i);
    t2 = ms_tbl.End(i);
    idxb = adv.time_avg >= t1 & adv.time_avg <= t2; % adv bursts
    idxy = edges.time >= t1 & edges.time <= t2; % adv edges
    % trim data to segment
    nb = sum(idxb);
    ti = adv.time_avg(idxb);
    u_adv = adv.vel_avg(idxb,:);
    y_advi = (edges.boundary_location(idxb)-157)/10;
    % calculate outer velocity from adcp
    u_outeri = nan(nb,4);
    for j = 1:nb % loop through bursts
        idxj = adcp.burst.time >= ti(j) & adcp.burst.time <= ti(j)+seconds(b_dt);
        velj = squeeze(mean(adcp.burst.vel_ice(idxj,idxr,[1 3 2]),2,'omitnan'));
        velj = [velj vecnorm(velj,2,2)];
        u_outeri(j,:) = mean(velj,'omitnan');
    end
    
    too_slow = [abs(u_outeri(:,1:3))<0.01 abs(u_outeri(:,4))<0.03];
    
    u_normi = u_adv./u_outeri;
    u_normi(too_slow) = nan;

    w_normi = u_adv(:,3)./abs(u_outeri(:,3));
    w_normi(too_slow(:,3)) = nan;

    % output to arrays
    y_adv = [y_adv; y_advi];
    u_norm = [u_norm; u_normi];
    w_norm = [w_norm; w_normi];
    u_outer = [u_outer; u_outeri];
end

idx_ice = y_adv < -1;
y_adv(idx_ice) = [];
u_norm(idx_ice,:) = [];
w_norm(idx_ice,:) = [];
u_outer(idx_ice,:) = [];

%% plot
figure(1); clf
clear ax
comp_lbls = {'u','v','w','|u|'};
for j = 4:-1:1
    ax(j) = subplot(4,1,j); hold on
    ylabel(ax(j),sprintf('%s_{adv}/%s_{adcp}',comp_lbls{j},comp_lbls{j}))
    xlabel(ax(j),'y [cm]')
end
figure(2); clf
ax_w = axes(figure(2)); hold on
ylabel(ax_w,'w_{adv}/|w|_{adcp}')
xlabel(ax_w,'y [cm]')
for i = 1:nsegs
    for j = 1:4
        plot(ax(j),y_adv,u_norm(:,j),'ko','markersize',3,'markerfacecolor','k')
        % errorbar(ax(j),y_adv,u_norm(:,j,1),u_norm(:,j,2)/10,'ko','linewidth',0.8,'markersize',3,'markerfacecolor','k')
    end
    plot(ax_w,y_adv,w_norm,'ko','markersize',3,'markerfacecolor','k')
    y_max = max([y_max max(y_adv)]);

end

linkaxes(ax,'x')
xlim(ax(1),[-0.5 max(ceil(2*y_max))/2])
xlim(ax_w,[-0.5 max(ceil(2*y_max))/2])