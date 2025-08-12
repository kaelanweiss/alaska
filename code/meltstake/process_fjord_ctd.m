% Script to collect CTD profiles from surface vessels (i.e. Polly and
% RiffRaft) during the 2024 meltstake deployments. CTD profiles are saved
% in a file in the processed data directory.
%
% KJW
% 4 Mar 2025

clear

% processed, combined ctd files
fname_p = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Polly\deployments_all\ctd_combo.mat';
fname_r = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\RiffRaft\combo_files\ctd_combo_RiffRaft.mat';

% processed meltstake directory
proc_dir = 'F:meltstake\data\proc';

% meltstake table
ms_tbl = loadMSInfo(28);

ctd_p = load(fname_p);
ctd_r = load(fname_r);

ctd_p = ctd_p.ctd;
ctd_r = ctd_r.ctd;

%% concatenate polly casts into a single structure like ctd_r
ctd_p_all = struct();
ndeps = length(ctd_p);
nc = 0;
for i = 1:ndeps
    nc = nc + length(ctd_p(i).time);
end
% preallocate fields
nz = length(ctd_p(1).depth);
flds = fieldnames(ctd_p);
flds = setxor(flds,{'depth','rsk'});
ctd_p_all.depth = ctd_p(1).depth;
for i = 1:length(flds)
    nr = size(ctd_p(1).(flds{i}),1);
    ctd_p_all.(flds{i}) = nan(nr,nc);
end
% fill in data
pos = 1;
for i = 1:ndeps
    % skip if empty
    if isempty(ctd_p(i).depth)
        continue
    end
    assert(all(size(ctd_p(i).depth)==size(ctd_p_all.depth)))
    % number of casts
    nci = length(ctd_p(i).time);
    
    % loop through fields
    for j = 1:length(flds)
        ctd_p_all.(flds{j})(:,pos:pos+nci-1) = ctd_p(i).(flds{j});
    end
    
    pos = pos + nci;
end

ctd_p = ctd_p_all;

% convert to datetime
ctd_p.time = dn2dt(ctd_p.time);
ctd_r.time = dn2dt(ctd_r.time);
ctd_p.time_full = dn2dt(ctd_p.time_full);
ctd_r.time_full = dn2dt(ctd_r.time_full);

%% combine ctd records
% set fields to extract
flds_extract = {'depth','time','T','S','PD','C','lon','lat','time_full'};

ctd_all(2) = struct();
ctd_all(1).vessel = 'polly';
ctd_all(2).vessel = 'riffraft';
for i = 1:length(flds_extract)
    fldi = flds_extract{i};
    ctd_all(1).(fldi) = ctd_p.(fldi);
    ctd_all(2).(fldi) = ctd_r.(fldi);
end

%% extract profiles during each meltstake deployment
deps_2024 = find(year(ms_tbl.Start)==2024)';
min_profs = 3; % minimum number of profiles to accept t_pad
for i = deps_2024
    fprintf('%d:',i)
    % default 30 min pad
    t_pad = hours(8.5);
    n_profs = 0;
    while n_profs < min_profs
        % time window
        t_pad = t_pad + hours(1);
        t1 = ms_tbl.Start(i)-t_pad;
        t2 = ms_tbl.End(i)+t_pad;
        % new copy of ctd_all to trim down and save
        ctd = ctd_all;
        % loop through vessels
        for j = 1:2
            idxt = ctd_all(j).time>=t1 & ctd_all(j).time<=t2;
%             fprintf(' %d',sum(idxt))
            % loop through fields
            for k = 2:length(flds_extract)
                fldk = flds_extract{k};
                ctd(j).(fldk) = ctd(j).(fldk)(:,idxt);
            end
            ctd(j).t_pad = t_pad;
        end
        % recalculate total number of profiles in window
        n_profs = length(ctd(1).time)+length(ctd(2).time);
    end

%     fprintf(' (%d hr)\n',hours(t_pad))
    
    % save file
    proc_fname = fullfile(proc_dir,ms_tbl.Folder{i},'ctd.mat');
    if ~exist(fullfile(proc_dir,ms_tbl.Folder{i}),'dir')
        mkdir(fullfile(proc_dir,ms_tbl.Folder{i}))
    end
    save(proc_fname,'ctd')
end




