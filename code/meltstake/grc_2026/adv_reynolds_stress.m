% Script to calculate Reynolds stresses from a full Meltstake deployment.
% This is a first crack at looking at near-boundary momentum fluxes.
%
%
% KJW
% 2 Jun 2026
clear

proc_dir = 'F:/meltstake/data/proc';
raw_dir = 'F:/meltstake/data/raw';

% deployment number
dep_num = 28;
ms_tbl = loadMSInfo(dep_num,'segments');
nsegs = size(ms_tbl,1);
dep_name = ms_tbl.Folder{1};

% load data
load(fullfile(proc_dir,dep_name,'adv','pck_edges.mat'))
load(fullfile(proc_dir,dep_name,'adv','svol_in_ocean.mat'))
load(fullfile(raw_dir,dep_name,'adv','adv.mat'),'adv')
load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'))
int_vel = load(fullfile(proc_dir,dep_name,'velocity_lowpass.mat'));

%% transform velocities (may be redundant)
switch dep_num
    case 26
        adv_roll_offset = 0;
    case 27
        adv_roll_offset = -21;
    case 28
        adv_roll_offset = 12.5;
end
adcp = msADCPTransform(adcp,adcp.burst.processing.cor_min,adcp.burst.processing.amp_min);
attitude = adcp.attitude;
attitude.roll_corrected = attitude.roll_corrected - adv_roll_offset;
adv = msADVTransform(adv,attitude);

% get rid of velocities in the ice
adv.vel_ice(~idx_ocean,:) = nan;

% nan out velocities not included in segments
idx_segs = false(size(adv.time));
for i = 1:nsegs
    t1 = ms_tbl.Start(i);
    t2 = ms_tbl.End(i);
    idx_segs = idx_segs | (adv.time>=t1 & adv.time<=t2);
end
adv.vel_ice(~idx_segs,:) = nan;

%% wall-normal coordinate adjustments based on flow direction
vel_seg_mean = nan(nsegs,6);
for i = 1:nsegs
    t1 = ms_tbl.Start(i);
    t2 = ms_tbl.End(i);
    % adcp
    idxt = adcp.burst.time>=t1 & adcp.burst.time<=t2;
    vel_seg_mean(i,1:3) = squeeze(mean(mean(adcp.burst.vel_ice(idxt,1:4,[1 3 2]),1,'omitnan'),2));
    % adv
    idxt = adv.time>=t1 & adv.time<=t2;
    vel_seg_mean(i,4:6) = mean(adv.vel_ice(idxt,:),'omitnan');
end
phi_adcp = atan2(vel_seg_mean(:,2),vel_seg_mean(:,1));
phi_adv = atan2(vel_seg_mean(:,5),vel_seg_mean(:,4));

% fudge 5 deg rotation for near-wall deployments, shouldn't affect further
% ones much
% adv_uv = adv.vel_ice(:,1)+1i*adv.vel(:,2);
% adv_uv = adv_uv/exp(1i*(-0)*pi/180);
% adv.vel_ice(:,1) = real(adv_uv);
% adv.vel_ice(:,2) = imag(adv_uv);

%% interpolate outer velocity onto adv time axis
tau_idx = 4;
tau = int_vel.tau(tau_idx); fprintf('Averaging timescale: %.1f s\n',tau)
adv.vel_outer = interp1(int_vel.time,int_vel.vel_lowpass(:,1,tau_idx+1)+1i*int_vel.vel_lowpass(:,3,tau_idx+1),adv.time);
adv.vel_angle = angle(adv.vel_outer);

% select vel range
vel_range = [.04 .06];
idx_vel_range = abs(adv.vel_outer)>vel_range(1) & abs(adv.vel_outer)<vel_range(2);
adv.vel_ice(~idx_vel_range,:) = nan;

% unwind velocity
adv_uw = adv.vel_ice(:,1) + 1i*adv.vel_ice(:,3);
adv_rot = adv_uw./exp(1i*adv.vel_angle);
adv.vel_rot = [real(adv_rot) adv.vel_ice(:,2) imag(adv_rot)];

%% Reynolds decomposition
% get ADV data on a continuous time vector
t = (adv.time(1):seconds(1/adv.sample_rate):adv.time(end))';
vel = nan(size(t,1),3);
burst_nums = unique(adv.burst_num);
orig_data = false(size(t));
for i = 1:length(burst_nums)
    idxi = adv.burst_num == burst_nums(i);
    ns = sum(idxi);
    t0i = adv.time(find(idxi,1));
    idx0 = find(t==t0i);
    if length(idx0) ~= 1
        warning('iter %d has %d starting indices',i,length(idx0))
    end
    vel(idx0:idx0+ns-1,:) = adv.vel_rot(idxi,:);
    orig_data(idx0:idx0+ns-1) = true;
end

% just a tiny bit of QC
vel(abs(vel)>0.8) = nan;

% smooth at a timescale and decompose
vel_bar = vel;
tau_rd = tau;
for i = 1:3
    vel_bar(:,i) = hannFilter(vel(:,i),adv.sample_rate*tau_rd);
end
vel_prm = vel(orig_data,:) - vel_bar(orig_data,:);

% correlations
corr = nan(size(vel_prm,1),3,3);
for i = 1:3
    for j = 1:3
        corr(:,i,j) = vel_prm(:,i).*vel_prm(:,j);
    end
end
adv.corr = [corr(:,2,1) corr(:,3,1) corr(:,3,2) corr(:,1,1) corr(:,2,2) corr(:,3,3)];
adv.tke = corr(:,1,1)+corr(:,2,2)+corr(:,3,3);

% average by burst
adv.vel_all_flds = {'outer_adcp','u_rot','v_rot','w_rot','uv','uw','vw','uu','vv','ww','tke'};
adv.vel_all = [adv.vel_outer adv.vel_rot adv.corr adv.tke];
adv = burstAverageADV(adv,'vel_all');
adv.y_avg = [nan; (edges.boundary_location-157)/10]; % this is not a good line of code

% normalized values
u_norm = adv.vel_avg(:,2:4)./abs(adv.vel_avg(:,1));
uv_norm = adv.vel_avg(:,5:7)./adv.vel_avg(:,11);
% uv_norm = adv.vel_avg(:,5:7)./abs(adv.vel_avg(:,1)).^2;
Cd = -adv.vel_avg(:,5)./abs(adv.vel_avg(:,1)).^2;
Cd(abs(adv.vel_avg(:,1))<.0 & adv.y_avg>500) = nan;
Cd(abs(Cd)>.1) = nan;

% bin in distance
[y_adv,sidx] = sort(adv.y_avg);
sidx(isnan(y_adv)) = [];
y_adv(isnan(y_adv)) = [];
[u_bin,y_bin] = time_bin(y_adv,[adv.vel_avg(sidx,1) u_norm(sidx,:) uv_norm(sidx,:) adv.vel_avg(sidx,8:11)]',0.5,'omitnan');
u_bin = u_bin';

% trim to y>0.55
idx_close = y_bin < 0.55;
u_bin(idx_close,:) = nan;

%% plots
fs = 12;
% velocity components
vel_lbls = {'along-flow','wall-normal','cross-flow'};
figure(11); clf
clear ax
for i = 3:-1:1
    ax(i) = subplot(3,1,i); hold on; box on
    plot(y_bin,u_bin(:,i+1),'ko','markersize',6,'markerfacecolor',colors(2))
    plot(adv.y_avg,u_norm(:,i),'kx')
    xlabel('y [cm]','fontsize',fs)
    ylabel({vel_lbls{i},'|u|_{ADCP}-normalized'},'fontsize',fs)
end
linkaxes(ax,'x')
xlim(ax(1),[0 12])
ylim(ax(1),[-.5 2])
ylim(ax(2),0.5*[-1 1])
ylim(ax(3),0.5*[-1 1])
title(ax(1),sprintf('%s | T_{outer}: %.1fs | T_{Rnlds}: %.1fs',strrep(dep_name,'_','\_'),tau,tau_rd),'fontsize',fs)

% reynolds stress components
vel_lbls = {'-<u''v''>/<q>','-<u''w''>/<q>','-<v''w''>/<q>'};
figure(12); clf
clear ax
for i = 3:-1:1
    ax(i) = subplot(3,1,i); hold on; box on
    plot(y_bin,-u_bin(:,i+4),'ko','markersize',6,'markerfacecolor',colors(2))
    plot(adv.y_avg,-uv_norm(:,i),'kx')
    xlabel('y [cm]','fontsize',fs)
    ylabel(vel_lbls{i},'fontsize',fs)
end
title(ax(1),sprintf('%s | T_{outer}: %.1fs | T_{Rnlds}: %.1fs',strrep(dep_name,'_','\_'),tau,tau_rd),'fontsize',fs)

% TKE
vel_lbls = {'<q>','<u''^2>','<v''^2>','<w''^2>'};
figure(13); clf
clear ax
for i = 4:-1:2
    ax(i) = subplot(4,1,i); hold on; box on
    plot(y_bin,u_bin(:,i+6),'ko','markersize',6,'markerfacecolor',colors(2))
    plot(adv.y_avg,adv.vel_avg(:,i+6),'kx')
    xlabel('y [cm]','fontsize',fs)
    ylabel(vel_lbls{i},'fontsize',fs)
end
ax(1) = subplot(4,1,1); hold on; box on
plot(adv.y_avg,adv.vel_avg(:,11),'kx')
plot(y_bin,u_bin(:,11),'ko','markersize',6,'markerfacecolor',colors(2))
xlabel('y [cm]','fontsize',fs)
ylabel(vel_lbls{1},'fontsize',fs)
title(ax(1),sprintf('%s | T_{outer}: %.1fs | T_{Rnlds}: %.1fs',strrep(dep_name,'_','\_'),tau,tau_rd),'fontsize',fs)
