function [T_scale,T_std,t_scale] = msTempScale(T,time,n,L)
% Function to calculate T scale to evaluate melt forcing based on RBR Solo
% records.
%
% [T_scale,T_std,t_scale] = msTempScale(T,time,n,L)
%
% Input
%   T: (nt x ns) matrix of temperature records
%   time: time vector corresponding to T
%   n: number of windows (use getWindowLength function)
%   L: length of windows in datetime units (use getWindowLength function)
%
% Output
%   
%
% KJW
% 30 Jan 2024
%
% NOT USED FOR ACTUAL T SCALE CALCULATION, ACTUAL SCRIPT IS
% "ms_scales_osm.m" - 30 Sep 2024

% convert datenum to datetime
if ~isdatetime(time)
    time = datetime(time,'convertfrom','datenum');
end

% values
ns = size(T,2);
t0 = time(1);

% preallocate
T_scale = nan(n,ns+1);
T_std = nan(n,ns+1);
t_scale = NaT(n,1);

% main calculation
for i = 1:n
    % time window
    t1 = t0 + (L/2)*(i-1);
    t2 = t1 + L;
    idx = time >= t1 & time < t2;
    t_scale(i) = t1 + L/2;

    % values
    % each thermistor
    for j = 1:ns
        Tj = T(idx,j);
        T_scale(i,j) = mean(Tj,'omitnan');
        T_std(i,j) = std(Tj,'omitnan');
    end
    
    % all thermistors
    Ti = T(idx,:);
    Ti = Ti(:);
    T_scale(i,ns+1) = mean(Ti,'omitnan');
    T_std(i,ns+1) = std(Ti,'omitnan');
end

