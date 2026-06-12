% Script to see if negative T is correlated with vertical velocity
% perturbations in the ADV data.
%
% KJW
% 28 May 2026
clear

proc_dir = 'F:/meltstake/data/proc';
raw_dir = 'F:/meltstake/data/raw';

% deployment number
dep_num = 28;
ms_tbl = loadMSInfo(dep_num,'segments');
nsegs = size(ms_tbl,1);
dep_name = ms_tbl.Folder{1};

% load
load(fullfile(raw_dir,dep_name,'rbr','T.mat'))
load(fullfile(raw_dir,dep_name,'adv','adv.mat'))
load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'))
load(fullfile(proc_dir,dep_name,'adv','svol_in_ocean.mat'))
load(fullfile(proc_dir,dep_name,'velocity_lowpass.mat'))
solo = T;

% transform
adv = msADVTransform(adv,adcp.attitude);

% trim
t1 = ms_tbl.Start(1);
t2 = ms_tbl.End(end);
idx_T1 = solo(1).time>=t1 & solo(1).time<=t2;
idx_T3 = solo(3).time>=t1 & solo(3).time<=t2;
idx_adv = adv.time >=t1 & adv.time <= t2;

t_T = solo(1).time(idx_T1);
t_adv = adv.time(idx_adv);
T = [solo(1).values(idx_T1) solo(3).values(idx_T3)];
dT = T(:,1)-T(:,2);
vel = adv.vel_ice;
vel(~idx_ocean,:) = nan;
vel = vel(idx_adv,:);

% smooth velocity onto T time grid, downsample
vel_ds = vel;
dt_solo = mean(seconds(diff(t_T)));
for i = 1:3
    vel_ds(:,i) = hannFilter(vel_ds(:,i),adv.sample_rate*dt_solo);
end
vel_ds = interp1(t_adv,vel_ds,t_T);


figure(1); clf
for i = 1:3
    subplot(1,3,i)
    plot(T(:,1),vel_ds(:,i),'k.')
end