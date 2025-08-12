% Script to create a processed CTD structure for each meltstake deployment 
% from raw RSK files.
%
% KJW
% 2 Feb 2024

clear

% data locations
raw_path = 'F:meltstake\data\raw\ctd';
proc_path = 'F:meltstake\data\proc';

% load meltstake deployment info
ms_tbl = loadMSInfo;
ndeps = size(ms_tbl,1);

% get list of available ctd files
d = dir(fullfile(raw_path,'*.rsk'));
rsk_files = cell(length(d),1);
for i = 1:length(d)
    rsk_files{i} = fullfile(d(i).folder,d(i).name);
end

% loop through deployments
tpad = hours(24);
dz = 0.5;
for i = 15%1:ndeps
    t1 = ms_tbl.Start(i) - tpad;
    t2 = ms_tbl.End(i) + tpad;
    
    fprintf('===== Deployment %d/%d: %s =====\n',i,ndeps,ms_tbl.Folder{i})
    try
        ctd = parseCTDProfiles(rsk_files,dz,[t1 t2]);
        %save(fullfile(proc_path,ms_tbl.Folder{i},'ctd.mat'),'ctd')
    catch
        fprintf('FAILED: %d/%d\n',i,ndeps)
    end
end
