% Script to work through what info to include in the TS plots of the
% boundary layer characteristics figure. Also figure out how to get rid of
% unphysical TS trajectories from the sensor misbehaving.
%
% KJW
% 24 Nov 2025
clear

raw_dir = 'F:/meltstake/data/raw';
proc_dir = 'F:/meltstake/data/proc';

load platform_comparison_data\ms_psd.mat

dep_tbl = loadMSInfo(26:28);
seg_tbl = loadMSInfo(26:28,'segments');

%% load hobo, T, polly ctd data
clear hobo TS
S0 = [26.7 27.8 28];
T0 = [6 6.9 7.1];
ms_depth = [21 43 52]; % m
TS(3) = struct();
for i = 1:3
    % hobo
    hobo_load = load(fullfile(raw_dir,dep_tbl.Folder{i},'hobo','hobo.mat'));
    hobo(i) = hobo_load.hobo(2);

    % TS
    idxt = hobo(i).time >= dep_tbl.Start(i) & hobo(i).time <= dep_tbl.End(i);
    TS(i).t = hobo(i).time(idxt);
    TS(i).S = hannFilter(hobo(i).S(idxt),1);
    TS(i).T = hannFilter(hobo(i).T_cal(idxt),1);
%     TS(i).fm = (T0(3)-Ti)/(T0(3)+90); % melt fraction
%     TS(i).fs = 1-Ti/T0(i); % fresh fraction
%     TS(i).rho0 = gsw_rho(S0(i),T0(i),ms_depth(i)); % outer density
%     TS(i).rho = gsw_rho(Si,Ti,ms_depth(i)); % bl density
%     TS(i).buoy = (9.8/rho0)*(rho0-rhoi); 
end

%% adcp data
% where to sample for velocities
y_s = [0.21] + 1e-6*[-1 1]; % distance from wall 
for i = 1:length(adcp)
    adcp(i).range = round(double(adcp(i).range),4);
    y_max = adcp(i).range(end);
    idx_s = (y_max-adcp(i).range) >= min(y_s) & (y_max-adcp(i).range) <= max(y_s);
    adcp(i).idx_s = idx_s;
    adcp(i).uvw = squeeze(mean(adcp(i).vel_ice(:,idx_s,:),2));
end
if ~(sum(adcp(1).idx_s) == sum(adcp(2).idx_s) && sum(adcp(1).idx_s) == sum(adcp(3).idx_s))
    warning('Number of bins do not agree')
end

%% clean up unphysical TS trajectories
k_hann = 3;
d_phi = 10*pi/180; % rad
fig = figure(100); clf
for i = 3
    Si = hannFilter(TS(i).S,k_hann);
    Ti = hannFilter(TS(i).T,k_hann);
    plot(Si,Ti,'k.-')
    hold on
    dT = diff(hannFilter(TS(i).T,k_hann));
    dS = diff(hannFilter(TS(i).S,k_hann));
    phi = atan2(dT,dS);
    phi = phi + pi/2;
    phi(phi<0) = phi(phi<0)+2*pi;
    phi = phi - pi/2;
    idx = abs(phi)<d_phi | abs(phi-pi)<d_phi;
    [blocks,widths] = findBlocks(idx);
    plot(Si(idx),Ti(idx),'r.')

end





