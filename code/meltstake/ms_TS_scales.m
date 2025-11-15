% Script to calculate temperature and salinity scales for meltstake
% segments. The script relies on the windows specified in the
% meltstake_deployment.xlsx spreadsheet as well as processed Solo and
% Concerto data. This script finally formalizes and centralizes the T and S
% scale methodology.
%
% KJW
% 18 Mar 2025
%
% IN PROGRESS...

clear

dep_nums = 28; % set to "nan" to load all
ms_tbl = loadMSInfo(dep_nums,'segments');

dep_names = unique(ms_tbl.Folder,'stable');
n_deps = length(dep_names);

% preallocate output vectors
n_all = size(ms_tbl,1);
depth = nan(n_all,1);
T = nan(n_all,2);
S = nan(n_all,2);
T_ctd = nan(n_all,2);
S_ctd = nan(n_all,2);
nS_ctd = nan(n_all,1);
pos = 0;

% loop through deployments
for i = 1:n_deps
    fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    
    dep_num = ms_tbl.Number(find(strcmp(ms_tbl.Folder,dep_names{i}),1));
    idx_dep = ms_tbl.Number == dep_num;
    proc_dir = fullfile('F:/meltstake/data/proc',dep_names{i});

    t_start = ms_tbl.Start(idx_dep);
    t_end = ms_tbl.End(idx_dep);

    % load ADCP data for depth
    load(fullfile('F:/meltstake/data/raw',dep_names{i},'adcp','adcp.mat'))
    adcp_depth = adcp.attitude.pressure-1;

    % load CTD data
    if exist(fullfile(proc_dir,'ctd.mat'),'file')
        load(fullfile(proc_dir,'ctd.mat'))
        includeS = 1;

        if length(ctd)==1 % pre-2024
            t_ctd = ctd.time;
            z_ctd = ctd.z;
            T_all_ctd = ctd.T;
            S_all_ctd = ctd.S;
        else % post-2024
            t_ctd = [ctd(1).time ctd(2).time];
            z_ctd = ctd(1).depth;
            T_all_ctd = [ctd(1).T ctd(2).T];
            S_all_ctd = [ctd(1).S ctd(2).S];
        end
    end

    % loop through windows in the deployment
    n_segs = sum(idx_dep);
    for j = 1:n_segs
        fprintf('segment %d\n',j)
        pos = pos + 1;
        t1 = t_start(j);
        t2 = t_end(j);
        
        % window depth
        depth(pos) = round(mean(adcp_depth(adcp.attitude.time>=t1 & adcp.attitude.time<=t2)),1);
        
        % temperature (onboard)
        T_file = fullfile(proc_dir,sprintf('T%d.mat',j));
        if exist(T_file,'file')
            T_solo = load(T_file);
            T(pos,:) = [mean(T_solo.T.T,'all','omitnan') std(T_solo.T.T(:),'omitnan')];
        end

        % salinity (onboard)
        S_file = fullfile(proc_dir,sprintf('hobo%d.mat',j));
        if exist(S_file,'file')
            S_hobo = load(S_file);
            S(pos,:) = [mean(S_hobo.hobo.S(:,2),'all','omitnan') std(S_hobo.hobo.S(:,2),'omitnan')];
        end

        % CTD salinity and temperature
        t_pad = 0.5*(t2-t1);
        if includeS
            idxt = t_ctd>=(t1-t_pad) & t_ctd<=(t2+t_pad);
            idxz = abs(z_ctd-round(2*depth(pos))/2)<=1;

            % expand search until CTD casts are found
            while sum(idxt)==0 && ( (t1-t_pad)>min(t_ctd) || (t2+t_pad)<max(t_ctd) )
                % find  preceding/following profiles
                t_pad = t_pad + minutes(10);
                idxt = t_ctd>=(t1-t_pad) & t_ctd<=(t2+t_pad);
                fprintf('expanding search window\n')
            end
            
            Sctdj = S_all_ctd(idxz,idxt);
            Tctdj = T_all_ctd(idxz,idxt);
            nS_ctd(pos) = numel(Sctdj);
            S_ctd(pos,:) = [mean(Sctdj,'all','omitnan') std(Sctdj(:),'omitnan')];
            T_ctd(pos,:) = [mean(Tctdj,'all','omitnan') std(Tctdj(:),'omitnan')];
            
        end
    end
end

%% print for entering in table
fprintf('depth:\n')
printScale(depth,1)

fprintf('\nT scale:\n')
printScale(T,1)

fprintf('\nS scale:\n')
printScale(S,1)

fprintf('\nT_ctd scale:\n')
printScale(T_ctd,1)

fprintf('\nS_ctd scale:\n')
printScale(S_ctd,1)

fprintf('\nnumber of S_ctd samples:\n')
printScale(nS_ctd,0)

function printScale(X,dec)
    n = size(X,1);
    m = size(X,2);
    fmt = sprintf('%%.%df (%%.%df)\n',dec,dec);
    for i = 1:n
        if m == 1
            fprintf([fmt(1:4) '\n'],round(X(i),dec))
        elseif m == 2
            fprintf(fmt,round(X(i,1),dec),round(X(i,2),dec))
        end
    end
end