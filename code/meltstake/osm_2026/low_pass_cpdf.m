% Script to redo cumulative pdf with low-pass filtering instead.
%
% KJW
% 11 Feb 2026
clear

raw_dir = 'F:/meltstake/data/raw';

% load
seg_tbl = loadMSInfo(28,'segments');
seg_tbl = seg_tbl(seg_tbl.include,:);
load(fullfile(raw_dir,seg_tbl.Folder{1},'adcp','adcp.mat'))

% trim to time window and max range
idxt = adcp.burst.time >= seg_tbl.Start(1) & adcp.burst.time <= seg_tbl.End(end);
r_min = min(seg_tbl.rmax);
idxr = adcp.burst.range <= (r_min-2*adcp.burst.cellsize);
vel = squeeze(mean(adcp.burst.vel_ice(idxt,idxr,[1 3 2]),2,'omitnan'));
vel = [vel sqrt(vel(:,1).^2+vel(:,2).^2)];
time = adcp.burst.time(idxt);

% timescales
fs = adcp.burst.samplerate;
T = seconds(time(end)-time(1));
tau = 2.^-(fix(log2(1/T)):fix(log2(fs)));

% calculate low passes
vel_lp = nan([size(vel) length(tau)]);
for i = 1:length(tau)
    w_hann = max([round(tau(i)*fs/2.5) 1]);
    for j = 1:size(vel,2)
        vel_lp(:,j,i) = hannFilter(vel(:,j),w_hann);
    end
end

% calculate KE and momentum
KE = nan([length(tau) 2]);
mom = nan([length(tau) 2]);
for i = 1:length(tau)
    KE(i,1) = mean(0.5*sum(vel_lp(:,1:3,i).^2,2));
    KE(i,2) = mean(0.5*vel_lp(:,4,i).^2);

    mom(i,1) = mean(sqrt(sum(vel_lp(:,1:3,i).^2,2)));
    mom(i,2) = mean(vel_lp(:,4,i));
end

