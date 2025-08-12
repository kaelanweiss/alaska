% Script to process and save T data in each meltstake deployment window.
% This takes data from the raw directory, cuts it into chunks according to
% "meltstake_deployments.xlsx", and saves the data to the processed 
% directory.
%
% KJW
% 22 Mar 2024

clear

ms_tbl = loadMSInfo('manualwindows');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');
n_deps = length(dep_nums);

%% loop through deployments
for i = 17:19%1:n_deps

    fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    idx_dep = dep_nums(i)==ms_tbl.Number;
    
    % load data
    raw_fldr = fullfile('F:meltstake\data\raw',dep_names{i});
    % Temp
    if exist(fullfile(raw_fldr,'rbr','T.mat'),'file')
        load(fullfile(raw_fldr,'rbr','T.mat'))
        time_all = dn2dt(T(1).time);
        T_all = repmat(T(1).values,[1 length(T)]);
        for j = 2:length(T)
            T_all(:,j) = T(j).values;
        end
    else
        warning('No temperature data found')
        continue
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
        T = struct('time',time_all(idxt),'T',T_all(idxt,:));
    
        % save trimmed and processed data structure
        proc_file = fullfile('F:meltstake\data\proc',dep_names{i},sprintf('T%d.mat',j));
        fprintf('saving %s\n',proc_file)
        save(proc_file,'T')
        fprintf('%.1f (%.1f)\n',round(mean(T.T,'all','omitnan'),1))
    
    end
end







