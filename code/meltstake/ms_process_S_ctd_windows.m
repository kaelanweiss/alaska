% Script to process and save S data in each meltstake deployment window.
% This takes data from the processed ctd file in the processed directory, 
% and finds an appropriate salinity value in the time frame specified in
% "meltstake_deployments.xlsx", and saves the data to the processed 
% directory. The scripts "ms_process_ctd.mat" and
% "ms_process_vel_windows.mat" should be run before so that processed ctd 
% and adcp (depth) files exist.
%
%
% KJW
% 24 Mar 2024

clear

% extend time window of section
t_pad = minutes(5);

ms_tbl = loadMSInfo('segments');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');
n_deps = length(dep_nums);

%% loop through deployments
for i = 1:n_deps

%     fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    idx_dep = dep_nums(i)==ms_tbl.Number;
    
    % load data
    ctd_file = fullfile('F:meltstake\data\proc',dep_names{i},'ctd.mat');
    if exist(ctd_file,'file')
        load(fullfile('F:meltstake\data\proc',dep_names{i},'ctd.mat'));
    else
        warning('File does not exist (deployment %d):\n%s',dep_nums(i),ctd_file)
        continue
    end
    
    % window
    n_wndws = sum(ms_tbl.Number==dep_nums(i));
    % loop through windows
    for j = 1:n_wndws
    
%         fprintf('--- Window %d ---\n',j)

        row_num = find(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        t1 = ms_tbl.Start(row_num);
        t2 = ms_tbl.End(row_num);

        % load adcp data to find depth
        load(fullfile('F:meltstake\data\proc',dep_names{i},['adcp' num2str(j) '.mat']));
        adcp_depth = adcp.attitude.pressure-1;
        zj = mean(adcp_depth(adcp.attitude.time>=t1 & adcp.attitude.time<=t2));
        
        % find window data
        idxt = ctd.t0>=(t1-t_pad) & ctd.t0<=(t2+t_pad);
        idxz = abs(ctd.z-round(2*zj)/2)<=1; % +/- 1 m  to the nearest 0.5 m
        
        if ~any(idxt)
            % find preceding and/or following profile(s)
            idx_prev = find(ctd.t0 < t1,1,'last');
            idx_next = find(ctd.t0 > t2,1,'first');
            idxt = [idx_prev idx_next];

%             fprintf('backward: %s\n',t1-ctd.t0(idx_prev))
%             fprintf('forward: %s\n',ctd.t0(idx_next)-t2)

            if isempty(idxt)
                warning('No ctd profiles found (deployment %d, window %d)',dep_nums(i),j)
                continue
            end            
        end
        
        % trim data
        ctd.z = ctd.z(idxz);
        ctd.time = ctd.time(idxz,idxt);
        ctd.t0 = ctd.t0(idxt);
        ctd.T = ctd.T(idxz,idxt);
        ctd.C = ctd.C(idxz,idxt);
        ctd.S = ctd.S(idxz,idxt);
        ctd.sigma = ctd.sigma(idxz,idxt);
        ctd.soundspeed = ctd.soundspeed(idxz,idxt);

        fprintf('%d\n',numel(ctd.S))
    
        % save trimmed and processed data structure
%         proc_file = fullfile('F:meltstake\data\proc',dep_names{i},sprintf('ctd%d.mat',j));
%         fprintf('saving %s\n',proc_file)
%         save(proc_file,'ctd')

%         figure(10*dep_nums(i)+j); clf
%         plot(ctd.S,ctd.z,'.-')
%         axis ij
%         if length(idxt)<=2
%             title('data outside window')
%         end
%         
    end
end







