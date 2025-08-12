% Script to start exploring unsteadiness in meltstake data
%
% KJW
% 6 Apr 2025

clear

addpath('..')

load F:/meltstake/data/raw/ms02_20240715_2001/adcp/adcp.mat

adcp = msADCPTransform(adcp,adcp.burst.processing.cor_min);

rmax = 0.4;
idxr = adcp.burst.range<=rmax;

% T_F = |u|/(d|u|/dt)

% bulk estimate using u_max
% smooth U
vel = adcp.burst.vel_ice(:,:,1:2); % u,w
% smooth across 3 points in the profile
for k = 1:size(vel,1)
    for m = 1:size(vel,3)
        vel(k,:,m) = hannFilter(squeeze(vel(k,:,m)),3);
    end
end
% calculate magnitude
vel_mag = vecnorm(vel,2,3);
% profile max and mean
[vel_max,i_max] = max(vel_mag,[],2);

% u_max
u_max = hannFilter(vel_max,1);

% u_mean
u_mean = mean(vel_mag(:,idxr),2);

% 3-pt derivative
dudt = gradient(u_max,1/adcp.burst.samplerate);
dudt(abs(dudt)<0.0001) = nan;

% unsteadiness time scale
T_Fu = abs(u_max./dudt);

% mean from all cells
vel_beam = adcp.burst.vel(:,adcp.burst.range<=rmax,:);

[~,dubdt,~] = gradient(vel_beam,1/adcp.burst.samplerate);
dubdt(abs(dubdt)<0.0001) = nan;
T_Fb = abs(vel_beam./dubdt);
