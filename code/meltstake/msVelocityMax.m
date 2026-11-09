function [vel_max,i_vel_max] = msVelocityMax(adcp,varargin)
% Function to calculate the maximum boundary layer forcing velocity used
% throughout the meltstake analysis. Smoothing in both space and time are
% performed before finding the maximum ice-parallel velocity.
%
% vel_max = msVelocityMax(adcp)
% vel_max = msVelocityMax(adcp,r_max)
%
% Input
%   adcp: adcp data structure, must contain burst.vel_ice
%   r_max: (optional) maximum adcp range to include
%
% Output
%   vel_max: time-dependent maximum ice-parallel velocity
%   i_vel_max: index of maximum velocity
%
% KJW
% 9 Nov 2026

% parse inputs
p = inputParser;
addOptional(p,'r_max',-1,@isnumeric)
parse(p,varargin{:});
r_max = p.Results.r_max;


% try to find approx ice edge if none is specified
if r_max == -1
    vel_mag = sqrt(adcp.burst.vel_ice(:,:,1).^2 + adcp.burst.vel_ice(:,:,2).^2);
    vel_mag_prof = mean(vel_mag,'omitnan');
    [vmp_max,i_max] = max(vel_mag_prof);

    while i_max < size(vel_mag,2)
        i_max = i_max+1;
        if vel_mag_prof(i_max)<0.2*vmp_max
            break
        end
    end
    fprintf('ice edge detected: %.2f m',round(double(adcp.burst.range(i_max)),2))
    if i_max == size(vel_mag,2)
        fprintf(' (end of profile)')
    end
    fprintf('\n')
else
    i_max = find(adcp.burst.range-r_max < 0,1);    
end

% trim down velocity
vel = adcp.burst.vel_ice(:,1:i_max+1,1:2);

% smooth in space
[nt,nc,nb] = size(vel);
for i = 1:nt
    for j = 1:nb
        vel(i,:,j) = hannFilter(squeeze(vel(i,:,j)),3);
    end
end

% smooth in time
k_hann = 1*adcp.burst.samplerate;
for i = 1:nc
    for j = 1:nb
        vel(:,i,j) = hannFilter(squeeze(vel(:,i,j)),k_hann);
    end
end

% compute velocity magnitude maximum
vel_mag = vecnorm(vel,2,3);
[vel_max,i_vel_max] = max(vel_mag,[],2);