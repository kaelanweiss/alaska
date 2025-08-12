% Script to plot output from parsed meltstake files.
%
% KJW
% 2 Oct 2024

clear

% data folder
raw_path = 'F:meltstake\data\raw';
% raw_path = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\raw\meltstake';
% raw_path = 'S:data/raw/meltstake/';

% deployment name
% dep = 'ms02_20240712_2020';
dep = 'ms02_20240715_2001';
% dep = 'ms01_20240717_0050';

% deployment path
dep_path = fullfile(raw_path,dep);

%% load files
% T
load(fullfile(dep_path,'rbr','T.mat'))

% ADCP
load(fullfile(dep_path,'adcp','adcp.mat'))

% ADV
load(fullfile(dep_path,'adv','adv.mat'))

% CTD
load(fullfile(dep_path,'ctd','ctd.mat'))

% HOBO
try
    load(fullfile(dep_path,'hobo','hobo.mat'))
end

%% plots
% adcp
adcp.burst.vel_ice = adcp.burst.vel_ice(:,:,1:3);
figure(1); clf
ax_adcp = adcpQuickPlot(figure(1),adcp,'vel',max(abs(extrema(adcp.burst.vel))),NaT,[0 1],4);

% T
figure(3); clf
hold on
for i = 1:3
    plot(T(i).time,T(i).values)
end
box on

ylabel('temp [C]')
title('Solo T')

% adv distance
amp = adv.pck_amp;
figure(4); clf
clear ax
for i = 1:3
    ax(i) = axes(figure(4),'position',axgridpos(1,3,i,[0.02,0.05,0.1,0.12]));
    box on
%     pcolor(repmat(adv.pck_time,[1 size(adv.pck_dist,2)]),adv.pck_dist,adv.pck_amp(:,:,i))
    pcolor(adv.pck_time,adv.pck_dist(1,:),adv.pck_amp(:,:,i)')
    shading flat
    clim(extrema(adv.pck_amp))
    hold on
    if i == 1
        ylabel('dist (mm)')
    end
    caxis([110 200])
end
cbar = colorbar;
cbar.Position = [0.65+.25+.01 .12 .02 .76];

linkaxes(ax)

% adv
figure(5); clf
plot(adv.time(1:4:end),adv.vel(1:4:end,:),'.-')
title('ADV vel')
ylabel('beam vel [m/s]')

% HOBO C
try
    figure(6); clf
    hold on
    for i = 1:3
        plot(hobo(i).time,hobo(i).C)
    end
    
    ylabel(sprintf('conductivity [%s]',hobo(1).units{1}))
    title('HOBO C')
end
