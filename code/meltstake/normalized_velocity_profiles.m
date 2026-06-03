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

%% interpolate outer velocity onto adv time axis
tau_idx = 6;
tau = int_vel.tau(tau_idx); fprintf('Averaging timescale: %.1f s\n',tau)
adv.vel_outer = interp1(int_vel.time,int_vel.vel_lowpass(:,1,tau_idx+1)+1i*int_vel.vel_lowpass(:,3,tau_idx+1),adv.time);
adv.vel_angle = angle(adv.vel_outer);
idx_slow = abs(adv.vel_outer)>.18;
adv.vel_outer(idx_slow) = nan;
% perform normalization
adv.vel_ice(~idx_ocean,:) = nan;
adv_uw = adv.vel_ice(:,1) + 1i*adv.vel_ice(:,3);
adv_norm = adv_uw./adv.vel_outer;
adv.vel_norm = [real(adv_norm) adv.vel_ice(:,2)./abs(adv.vel_outer) imag(adv_norm) adv.vel_ice(:,4)./abs(adv.vel_outer)];
adv = burstAverageADV(adv,'vel_norm');

% normalized stresses


% figure to check that things are working
figure(11); clf
subplot(2,1,1); hold on
plot(adv.time,real(adv_uw))
plot(adv.time,real(adv.vel_outer))
ylim(0.4*[-1 1])
subplot(2,1,2); hold on
plot(adv.time,imag(adv_uw))
plot(adv.time,imag(adv.vel_outer))
ylim(0.4*[-1 1])

%% calculate velocity things
% set up arrays
y_adv = [];
u_norm = [];

y_max = 0;
for i = 1:nsegs
    fprintf('segment %d\n',i)
    % create indexing vectors
    t1 = ms_tbl.Start(i);
    t2 = ms_tbl.End(i);
    idxb = adv.time_avg >= t1 & adv.time_avg <= t2; % adv bursts
    idxy = edges.time >= t1 & edges.time <= t2; % adv edges
    idxt = adv.time >= t1 & adv.time <= t2; % adv full resolution
    % trim data to segment
    nb = sum(idxb);
    ti = adv.time_avg(idxb);
    u_normi = adv.vel_avg(idxb,:);
    y_advi = (edges.boundary_location(idxb)-157)/10;

    % output to arrays
    y_adv = [y_adv; y_advi];
    u_norm = [u_norm; u_normi];
end

idx_ice = y_adv < -1;
y_adv(idx_ice) = [];
u_norm(idx_ice,:) = [];

% sort
[y_adv,sidx] = sort(y_adv);
u_norm = u_norm(sidx,:);

% get rid of nans in y_adv and then bin by y
idx_nan = isnan(y_adv);
y_adv(idx_nan) = [];
u_norm(idx_nan,:) = [];
[u_bin,y_bin] = time_bin(y_adv,u_norm',0.5,'omitnan');
u_bin = u_bin';

%% plot
fs = 12;
figure(tau_idx+10); clf
clear ax
comp_lbls = {'along-flow','wall-normal','cross-flow','magnitude'};
for j = 4:-1:1
    ax(j) = subplot(4,1,j); hold on; box on
    ylabel(ax(j),sprintf('%s',comp_lbls{j}),'fontsize',fs)
    xlabel(ax(j),'y [cm]','fontsize',fs)
end

for i = 1:nsegs
    for j = 1:4
        plot(ax(j),y_adv,u_norm(:,j),'kx','markersize',5,'markerfacecolor','k')
        plot(ax(j),y_bin,u_bin(:,j),'ko','markersize',6,'markerfacecolor',colors(2))
    end
    y_max = max([y_max max(y_adv)]);
end
title(ax(1),sprintf('%s (T_{filter}=%.1f s)',strrep(dep_name,'_','\_'),tau),'fontsize',fs+1)
linkaxes(ax,'x')
xlim(ax(1),[-0.5 max(ceil(2*y_max))/2])

switch dep_num
    case 28
        ylim(ax(1),[-0.1 1.5])
        ylim(ax(2),[-0.3 0.3])
        ylim(ax(3),[-0.7 0.7])
        ylim(ax(4),[0 2])
end