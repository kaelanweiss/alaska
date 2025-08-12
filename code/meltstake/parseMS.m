% Script to parse meltstake data in the field. This takes files from the
% instruments and turns them into usable formats.
%
% KJW
% 16 Jun 2024

clear

% data folder
% raw_path = 'F:AK_202406/data/meltstake';
% raw_path = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\raw\meltstake';
% raw_path = 'S:data/raw/meltstake/';
raw_path = 'F:\meltstake\data\raw';

tbl_path = 'G:\Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','overview');
dep_num = 20;

% deployment name
% dep = 'ms01_20230708_0218';
dep = ms_tbl.Folder{dep_num};

% deployment path
dep_path = fullfile(raw_path,dep);

% set time range
% t1 = datetime(2023,7,8,1,46,0);
% t2 = datetime(2023,7,8,2,20,0);
t1 = ms_tbl.Start(dep_num);
t2 = ms_tbl.End(dep_num);

% save flag
save_files = 0;

%% parse data
% T
T = parseRBRDeployment(dep_path,t1,t2);

% ADV
[adv,raw] = parseADVDeployment(fullfile(dep_path,'adv'),t1,t2);

% CTD
ctd = parseRBRDeployment(dep_path,t1-minutes(15),t1+minutes(15),'ctd');
% ctd = parseRBRDeployment(dep_path,t1,t2,'ctd');

% HOBO
hobo = parseHOBODeployment(dep_path,t1-minutes(15),t2);

%% ADCP
cor_min = 50;
adcp = parseNortekDeployment(fullfile(dep_path,'adcp'),t1,t2);
adcp = msADCPTransform(adcp,cor_min,0);

adcpQuickPlot(figure(1),adcp,'vel',0.18*[-1 1],NaT,[0 1.2],8);
adcpQuickPlot(figure(2),adcp,'vel_ice',0.30*[-1 1],NaT,[0 1.2],8);

%% plots
% T
figure(3); clf
hold on
for i = 1:3
    plot(T(i).time,T(i).values)
end

ylabel('temp [C]')
title('Solo T')

% adv distance
amp = adv.pck_amp;
figure(4); clf
clear ax
for i = 1:3
    ax(i) = axes(figure(4),'position',axgridpos(1,3,i,[0.02,0.05,0.1,0.12]));
    box on
    pcolor(repmat(adv.pck_time,[1 size(adv.pck_dist,2)]),adv.pck_dist,adv.pck_amp(:,:,i))
    shading flat
    clim(extrema(adv.pck_amp))
    hold on
    if i == 1
        ylabel('dist (mm)')
    end
    caxis([100 200])
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
figure(6); clf
hold on
for i = 1:3
    plot(hobo(i).time,hobo(i).C)
end

ylabel(sprintf('conductivity [%s]',hobo(1).units{1}))
title('HOBO C')

%% save files
if save_files
    % T
    save(fullfile(dep_path,'rbr','T.mat'),'T')

    % ADCP
    save(fullfile(dep_path,'adcp','adcp.mat'),'adcp','-v7.3')

    % ADV
    save(fullfile(dep_path,'adv','adv.mat'),'adv','raw')

    % CTD
    save(fullfile(dep_path,'ctd','ctd.mat'),'ctd')

    % HOBO
    save(fullfile(dep_path,'hobo','hobo.mat'),'hobo')

end
