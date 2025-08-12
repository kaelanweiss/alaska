% Script to make some quick spectra of velocity at the glacier face.
%
% KJW
% 13 Jan 2025

clear

raw_path = 'F:AK_202406/data/meltstake';

dep_names = {'ms02_20240712_2020','ms02_20240715_2001','ms01_20240717_0050'};

% time windows in each deployment
tlims = [datetime(2024,7,12,20,43,0) datetime(2024,7,12,21,26,0);...
         datetime(2024,7,15,20,2,0) datetime(2024,7,15,22,58,0);...
         datetime(2024,7,17,1,12,0) datetime(2024,7,17,3,21,0)];

% max range
rmax = [0.5 0.4 0.5];

%% load
clear adcp
for i = 1:length(dep_names)
    load_temp = load(fullfile(raw_path,dep_names{i},'adcp','adcp.mat'));
    adcp(i) = load_temp.adcp;
end

%% psds
clear psd
psd(length(adcp)) = struct('f',[],'P',[]);

for i = 1:length(adcp)
    % sample frequency
    fs = adcp(i).burst.samplerate; % s^-1
    % mask data by time and range
    idxt = adcp(i).burst.time >= tlims(i,1) & adcp(i).burst.time <= tlims(i,2);
    idxr = adcp(i).burst.range <= rmax(i);
    % extract data
    vel = adcp(i).burst.vel(idxt,idxr,1:5);
    [nt,nr,nb] = size(vel);
    % patch holes (this is a bad way to do this)
    vel(isnan(vel)) = 0;
    % preallocate
    nf = ceil(nt/2)+1;
    psd(i).P = nan(nf,nr,nb);
    % loop through beams
    for j = 1:nb
        [Pj,f] = psd_matlab(vel(:,:,j),fs);
        psd(i).P(:,:,j) = Pj;
    end
    psd(i).f = f;
    psd(i).P_mean = squeeze(mean(psd(i).P,2));
end

%% plots
figure(1); clf





