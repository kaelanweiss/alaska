% KJW
% 9 Mar 2025

clear

% data folder
raw_path = 'F:meltstake\data\raw';
% raw_path = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\raw\meltstake';
% raw_path = 'S:data/raw/meltstake/';

% deployment name
% dep = 'ms02_20240712_2020';
dep = 'ms02_20240715_2001';
% dep = 'ms01_20240717_0050';

ms_tbl = loadMSInfo;
row_num = strcmp(ms_tbl.Folder,dep);
t1 = ms_tbl.Start(row_num);
t2 = ms_tbl.End(row_num);

% deployment path
dep_path = fullfile(raw_path,dep);

%% load files
% T
load(fullfile(dep_path,'rbr','T.mat'))

% ADCP
load(fullfile(dep_path,'adcp','adcp.mat'))
adcp.burst.vel_ice(:,:,1) = -adcp.burst.vel_ice(:,:,1);
adcp.burst.vel_ice(:,:,3) = -adcp.burst.vel_ice(:,:,3);

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
adcp.burst.vel_ice = adcp.burst.vel_ice(:,:,1:4);
figure(1); clf
ax = adcpQuickPlot(figure(1),adcp,'vel_ice',0.25,[t1 t2],[0 0.8],4);
linkaxes(ax,'off')

% T
cla(ax(4))
hold on
for i = 1:3
    plot(ax(4),T(i).time,T(i).values)
end
box on
ylabel(ax(4),'temp [\circC]','fontsize',11)
ylim([2.5 6])
legend({'near','mid','far'},'location','south','orientation','horizontal','fontsize',10)

% title('Solo T')


