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
% raw_path = 'F:\meltstake\data\raw';
raw_path = 'F:/GIANT2026/meltstake/data/raw';

% CAREFUL RUNNING THIS ON OLD DEPLOYMENTS, IT WILL DELETE UNWRAPPED DATA
% dep_num = 26;
% ms_tbl = loadMSInfo(dep_num);

% deployment name
% dep = '20260730_1530_ms01';
% dep = '20260803_1430_ms01';
dep = '20260816_1850_ms01';
% dep = ms_tbl.Folder{1};

% deployment path
dep_path = fullfile(raw_path,dep);

% set time range
% t1 = datetime(2000,1,1);
% t2 = datetime('now');
t1 = datetime(2026,8,16,18,55,0);
t2 = datetime(2026,8,17,14,15,0);
% t1 = ms_tbl.Start;
% t2 = ms_tbl.End;

% save flag
save_files = 0;
%% parse data
fprintf('Parsing data for deployment %s\n',dep)
fprintf('t1: %s\nt2: %s\n',t1,t2)
fprintf('Saving files: %d\n',save_files)
% T
fprintf('Solos...\n')
T = parseRBRDeployment(dep_path,t1,t2,'solo');

% MS CTD
fprintf('MS CTD...\n')
ctd_ms = parseMSCTDDeployment(fullfile(dep_path,'ctd_ms'),t1,t2);

% ADV
fprintf('ADV...\n')
adv = parseADVDeployment(fullfile(dep_path,'adv'),t1,t2);

% CTD
fprintf('ROV CTD...\n')
ctd_rov = parseRBRDeployment(dep_path,t1-minutes(20),t1+minutes(20),'ctd_rov');
% ctd_rov = parseRBRDeployment(dep_path,t1,t2,'ctd_rov');

% HOBO
% hobo = parseHOBODeployment(dep_path,t1-minutes(15),t2);

% 881a
fprintf('Sonar...\n')
sonar = parse881aDeployment(fullfile(dep_path),t1,t2);

%% ADCP
fprintf('ADCP...\n')
cor_min = 50;
adcp = parseNortekDeployment(fullfile(dep_path,'adcp'),t1,t2);
adcp = msADCPTransform(adcp,cor_min,0);

adcpQuickPlot(figure(1),adcp,'vel',0.05*[-1 1],NaT,[0 0.8],8);
adcpQuickPlot(figure(2),adcp,'vel_ice',0.08*[-1 1],NaT,[0 0.8],8);
adcpQuickPlot(figure(30),adcp,'amp',[84 86],NaT,[0 1.2],8);

%% ADV transform
adv = msADVTransform(adv,adcp.attitude);

%% plots
% T
figure(3); clf
hold on
if exist('ctd_ms','var')
    plot(ctd_ms.time,ctd_ms.T,'k-')
end
for i = 1:3
    plot(T(i).time,T(i).values,'color',colors(i))
end
ylabel('temp [C]')
title('Solo/CTD T')

% S
figure(4); clf
plot(ctd_ms.time,ctd_ms.S)
ylabel('CTD salinity [PSU]')

% adv distance
figure(5); clf
clear ax
dt = 10;
for i = 1:3
    ax(i) = axes(figure(5),'position',axgridpos(1,3,i,[0.02,0.05,0.1,0.12]));
    box on
    pcolor(adv.pck_time(1:dt:end),adv.pck_dist(1,:),adv.pck_amp(1:dt:end,:,i)')
    shading flat
    clim(extrema(adv.pck_amp))
    hold on
    if i == 1
        ylabel('dist (mm)')
    end
    clim([100 200])
end
cbar = colorbar;
cbar.Position = [0.65+.25+.01 .12 .02 .76];

linkaxes(ax)

% adv
figure(6); clf
plot(adv.time(1:10:end),adv.vel(1:10:end,:),'.-')
title('ADV vel')
ylabel('beam vel [m/s]')

% sonar881a
n_scan = floor(length(sonar)/2);
ax_881a = plot881aScan(figure(7),sonar(n_scan));


%% save files
if save_files
    fprintf('Saving files...\n')
    % T
    if exist('T','var')
        fprintf('\t%s\n','solo/T.mat')
        save(fullfile(dep_path,'solo','T.mat'),'T')
    end

    % MS CTD
    if exist('ctd_ms','var')
        fprintf('\t%s\n','ctd_ms/ctd_ms.mat')
        save(fullfile(dep_path,'ctd_ms','ctd_ms.mat'),'ctd_ms')
    end

    % ADCP
    if exist('adcp','var')
        fprintf('\t%s\n','adcp/adcp.mat')
        save(fullfile(dep_path,'adcp','adcp.mat'),'adcp','-v7.3')
    end

    % ADV
    if exist('adv','var')
        fprintf('\t%s\n','adv/adv.mat')
        save(fullfile(dep_path,'adv','adv.mat'),'adv')
    end

    % ROV CTD
    if exist('ctd_rov','var')
        fprintf('\t%s\n','ctd_rov/ctd_rov.mat')
        save(fullfile(dep_path,'ctd_rov','ctd_rov.mat'),'ctd_rov')
    end

    % HOBO
    if exist('hobo','var')
        fprintf('\t%s\n','hobo/hobo.mat')
        save(fullfile(dep_path,'hobo','hobo.mat'),'hobo')
    end

    % 881a SONAR
    if exist('sonar','var')
        fprintf('\t%s\n','sonar881a/sonar881a.mat')
        save(fullfile(dep_path,'sonar881a','sonar881a.mat'),'sonar','-v7.3')
    end

end
