% Script to create an index vector for a meltstake deployment to determine
% when the ADV is sampling in water and not in ice. The intention is to
% use this somewhat manually built vector for further ADV processing.
%
% KJW
% 18 May 2026
clear

raw_dir = 'F:/meltstake/data/raw';
proc_dir = 'F:/meltstake/data/proc';

% one deployment at a time (will break otherwise)
ms_tbl = loadMSInfo(26,'segments');

% load adv data
dep_name = ms_tbl.Folder{1};
load(fullfile(raw_dir,dep_name,'adv','adv.mat'),'adv')

% set up indexing vector
idx_ocean = true(size(adv.time));

% set times in ice
switch ms_tbl.Number(1)
    case 26 % dep26
        % in_ice_list = [ms_tbl.Start(1) datetime(2024,7,12,20,34,6);...
        %           ms_tbl.Start(2) ms_tbl.End(2);...
        %           ms_tbl.Start(3) ms_tbl.End(3);...
        %           ms_tbl.Start(3) datetime(2024,7,12,21,34,30)];
        in_ice_list = {'12-Jul-2024 20:22:30' '12-Jul-2024 20:36:06';...
                       '12-Jul-2024 20:41:30' '12-Jul-2024 21:01:25';...
                       '12-Jul-2024 21:01:25' '12-Jul-2024 21:21:25';...
                       '12-Jul-2024 21:01:25' '12-Jul-2024 21:34:30'};
        
    case 27 % dep27
        % in_ice_list = [ms_tbl.Start(2) ms_tbl.End(2);...
        %           ms_tbl.Start(4) datetime(2024,7,15,21,08,05);...
        %           ms_tbl.Start(8) ms_tbl.End(8);...
        %           ms_tbl.Start(9) datetime(2024,7,15,22,35,05);...
        %           ];
        in_ice_list = {'15-Jul-2024 20:21:15' '15-Jul-2024 20:41:00';...
                       '15-Jul-2024 21:01:15' '15-Jul-2024 21:08:05';...
                       '15-Jul-2024 22:21:20' '15-Jul-2024 22:29:30';...
                       '15-Jul-2024 22:29:30' '15-Jul-2024 22:35:05'};

    case 28 % dep28
        % in_ice_list = [ms_tbl.Start(2) datetime(2024,7,17,01,24,00);...
        %           ms_tbl.Start(3) datetime(2024,7,17,01,45,15);...
        %           ];
        in_ice_list = {'17-Jul-2024 01:12:05' '17-Jul-2024 01:24:00';...
                       '17-Jul-2024 01:32:00' '17-Jul-2024 01:45:15'};

end

% loop through times above
for i = 1:size(in_ice_list,1)
    t1 = datetime(in_ice_list{i,1});
    t2 = datetime(in_ice_list{i,2});
    idxi = adv.time>=t1 & adv.time<=t2;
    idx_ocean(idxi) = false;
end

% save
save_output = 1;
if save_output
    save(fullfile(proc_dir,dep_name,'adv','svol_in_ocean.mat'),'idx_ocean','in_ice_list')
end