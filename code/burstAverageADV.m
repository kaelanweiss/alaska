function adv = burstAverageADV(adv,varargin)
% Function to calculate burst-averaged velocities in ADV data.
%
% adv = burstAverageADV(adv)
%
% Input
%   adv: ADV structure
%
% Output
%   adv: ADV structure with burst-averaged fields added
% 
% KJW
% 31 Jul 2025

% parse
fld = 'vel_ice';
if nargin > 1
    fld = varargin{1};
end

% get list of unique bursts
bursts = unique(adv.burst_num);
n_bursts = length(bursts);

% preallocate
nc = size(adv.(fld),2);
vel_avg = nan(n_bursts,nc);
vel_std = nan(n_bursts,nc);
n_avg = nan(n_bursts,1);
time_avg = NaT(n_bursts,1);

% loop through bursts
for i = 1:n_bursts
    % slice according to burst number
    bnum = bursts(i);
    idx = adv.burst_num == bnum;
    
    % get/calculate values
    time_avg(i) = adv.time(find(idx,1));
    n_avg(i) = sum(~isnan(adv.(fld)(idx,1)));
    vel_avg(i,:) = mean(adv.(fld)(idx,:),'omitnan');
    vel_std(i,:) = std(adv.(fld)(idx,:),'omitnan');
end

% add to ADV structure
adv.time_avg = time_avg;
adv.n_avg = n_avg;
adv.vel_avg = vel_avg;
adv.vel_std = vel_std;
adv.fld_avg = fld;
