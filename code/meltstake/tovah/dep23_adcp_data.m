% Script to clean up and export ADCP data from Spireberg deployment 23
%
% KJW
% 8 May 2026
clear

dep_tbl = loadMSInfo(23);
% load(fullfile('F:/meltstake/data/raw',dep_tbl.Folder{1},'adcp','adcp.mat'))

adcp = parseNortekDeployment(fullfile('F:/meltstake/data/raw',dep_tbl.Folder{1},'adcp'),dep_tbl.Start(1),dep_tbl.End(1));

cor_min = 30;
amp_min = 0;
adcp = msADCPTransform(adcp,cor_min,amp_min);

% smooth
vel = adcp.burst.vel_ice(:,:,[1 3 2]);
w_hann = 8; % sec
[nt,nc,nb] = size(vel);
for i = 1:nc
    for j = 1:nb
        vel(:,i,j) = hannFilter(vel(:,i,j),w_hann*adcp.burst.samplerate);
    end
end
adcp.burst.vel_smth = vel;

% volume-average
r_lim = [0.4 0.6];
idxr = adcp.burst.range >= r_lim(1) & adcp.burst.range <= r_lim(2);
vel_mean = squeeze(mean(vel(:,idxr,:),2,'omitnan'));

% setup data
t0 = adcp.burst.time(1);
t0_str = sprintf('%s',t0);
time = adcp.burst.time;
time_sec = seconds(adcp.burst.time-t0);
range = adcp.burst.range;
vel_uvw = vel;
vel_uvw_mean = vel_mean;

save dep23_adcp.mat t0_str time time_sec range vel_uvw vel_uvw_mean