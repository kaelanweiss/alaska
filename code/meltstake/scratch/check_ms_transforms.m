clear

raw_dir = 'F:meltstake/data/raw';
proc_dir = 'F:meltstake/data/proc';

% choose deployment and segment(s)
dep_num = 28;
seg_num = [1 1];
ms_tbl = loadMSInfo(28,'segments');

% find row number and time window
t1 = ms_tbl.Start(seg_num(1));
t2 = ms_tbl.End(seg_num(2));

% load
dep_name = ms_tbl.Folder{1};
load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'))
load(fullfile(raw_dir,dep_name,'adv','adv.mat'))
load(fullfile(proc_dir,dep_name,'adv','svol_in_ocean.mat'))

% trim data


% TEST %
pitch = 90;
roll = 180;
adcp.attitude.roll(:) = pitch+270;
adcp.attitude.pitch(:) = -roll;

% transform adcp velocity
% adcp.burst.vel_unw(:,:,5) = 0.1*adcp.burst.vel_unw(:,:,5);
adcp = msADCPTransform(adcp,adcp.burst.processing.cor_min,adcp.burst.processing.amp_min);

% time indexing
idxt_adcp = adcp.burst.time>=t1 & adcp.burst.time<=t2;
idxt_adv = adv.time>=t1 & adv.time<=t2;

% transform ADV velocity
adv = msADVTransform(adv,adcp.attitude);
adv.vel_ice(~idx_ocean,:) = nan;