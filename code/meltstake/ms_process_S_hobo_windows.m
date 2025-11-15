% Script to process and save hobo data in each meltstake deployment window.
% This takes data from the raw directory, cuts it into chunks according to
% "meltstake_deployments.xlsx", and saves the data to the processed 
% directory.
%
% KJW
% 10 Nov 2025

clear

ms_tbl = loadMSInfo('segments');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');
n_deps = length(dep_nums);

%% loop through deployments
for i = 1:n_deps

    fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    idx_dep = dep_nums(i)==ms_tbl.Number;
    
    % load data
    raw_fldr = fullfile('F:meltstake\data\raw',dep_names{i});
    % hobo file
    if exist(fullfile(raw_fldr,'hobo','hobo.mat'),'file')
        load(fullfile(raw_fldr,'hobo','hobo.mat'))
        hobo_all = hobo;
    else
        warning('No hobo data found')
        continue
    end

    % unified time axis
    time_all = hobo_all(1).time;
    nt = length(time_all);
    nh = length(hobo);
    T = nan(nt,nh);
    S = nan(nt,nh);
    for j = 1:nh
        T(:,j) = interp1(hobo(j).time,hobo(j).T_cal,time_all,'nearest');
        S(:,j) = interp1(hobo(j).time,hobo(j).S,time_all,'nearest');
    end
    
    % window
    n_wndws = sum(ms_tbl.Number==dep_nums(i));
    % loop through windows
    for j = 1:n_wndws
    
        fprintf('--- Window %d ---\n',j)
        
        row_num = find(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        t1 = ms_tbl.Start(row_num);
        t2 = ms_tbl.End(row_num);
        
        % find window data
        idxt = time_all>=t1 & time_all<=t2;
        
        % trim data down to time and range for each window
        hobo = struct('time',time_all(idxt),'T',T(idxt,:),'S',S(idxt,:));
    
        % save trimmed and processed data structure
        proc_file = fullfile('F:meltstake\data\proc',dep_names{i},sprintf('hobo%d.mat',j));
        fprintf('saving %s\n',proc_file)
        save(proc_file,'hobo')
    
    end
end







