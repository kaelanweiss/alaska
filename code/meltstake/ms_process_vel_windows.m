% Script to process and save ADCP data in each meltstake deployment window.
% This takes data from the raw directory, cuts it into chunks according and
% processes it with values specified in "meltstake_deployments.xlsx", and
% saves the data to the processed directory.
%
% KJW
% 22 Mar 2024

clear

ms_tbl = loadMSInfo(26,'manualwindows');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');
r_max = ms_tbl.rmax;
n_deps = length(dep_nums);

% times to exclude due to bad data, corresponding deployment number
t_exclude = [datetime(2023,5,29,21,59,0) datetime(2023,5,29,22,03,0);...
             datetime(2023,5,30,1,41,0) datetime(2023,5,30,1,47,0);...
             datetime(2023,7,9,21,30,0) datetime(2023,7,9,21,40,0);...
             datetime(2023,7,9,21,45,0) datetime(2023,7,9,21,51,30);...
             datetime(2023,7,9,23,22,0) datetime(2023,7,9,23,28,0);...
             datetime(2024,7,17,2,11,55) datetime(2024,7,17,2,11,58)];
t_exclude_dep = [2;...
                 3;...
                 8;...
                 8;...
                 8;...
                 28];

%% loop through deployments
for i = 1:n_deps

    % time window
    fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    idx_dep = dep_nums(i)==ms_tbl.Number;
    
    % load data
    load(fullfile('F:meltstake\data\raw',dep_names{i},'adcp','adcp.mat'))
    adcp_all = adcp;
    adcp_all.burst.time = dn2dt(adcp_all.burst.time);
    adcp_all.attitude.time = dn2dt(adcp_all.attitude.time);
    [nt,nc,nb] = size(adcp_all.burst.vel);
    
    % window
    n_wndws = sum(ms_tbl.Number==dep_nums(i));
    % loop through windows
    for j = 1:n_wndws
    
        fprintf('--- Window %d ---\n',j)
        
        row_num = find(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        t1 = ms_tbl.Start(row_num);
        t2 = ms_tbl.End(row_num);
        
        % process data
        cor_min = ms_tbl.cor_min(row_num);
        amp_min = 0;
        
        % find window data
        idxt = adcp_all.burst.time>=t1 & adcp_all.burst.time<=t2;
        idxr = adcp_all.burst.range<=r_max(row_num);
        
        % trim data down to time and range for each window
        adcp = adcp_all;
        if isfield(adcp,'echo')
            adcp = rmfield(adcp,'echo');
        end
        adcp.burst.time = adcp.burst.time(idxt);
        adcp.burst.range = adcp.burst.range(idxr);

        adcp.burst.vel = adcp.burst.vel(idxt,idxr,:);
        adcp.burst.cor = adcp.burst.cor(idxt,idxr,:);
        adcp.burst.amp = adcp.burst.amp(idxt,idxr,:);
        if isfield(adcp.burst,'vel_unw')
            adcp.burst.vel_unw = adcp.burst.vel_unw(idxt,idxr,:);
        end

        % exclude specified times with bad/unwanted data
        time = adcp.burst.time;
        idx_exclude = false(size(time));
        for k = 1:size(t_exclude,1)
            if t_exclude_dep(k)==dep_nums(i)
                idxk = time>=t_exclude(k,1) & time<=t_exclude(k,2);
                idx_exclude = idx_exclude & ~idxk;
            end
        end
        adcp.burst.vel(idx_exclude,:,:) = nan;
        if isfield(adcp.burst,'vel_unw')
            adcp.burst.vel_unw(idx_exclude,:,:) = nan;
        end

        % velocity transform
        adcp = msADCPTransform(adcp,cor_min,amp_min);

        % calculate u_max and u_coarse
        vel = adcp.burst.vel_ice(:,:,1:2); % u,w
        % smooth across 3 points in the profile
        for k = 1:size(vel,1)
            for m = 1:size(vel,3)
                vel(k,:,m) = hannFilter(squeeze(vel(k,:,m)),3);
            end
        end
        % calculate magnitude
        vel_mag = vecnorm(vel,2,3);
        % profile max
        adcp.burst.u_max = max(vel_mag,[],2);
        % coarse velocity
        uw_coarse = mean(squeeze(vel(:,2,:)),'omitnan');
        u_coarse = vecnorm(uw_coarse,2,2);
        uw_std = [std(vel(:,2,1),'omitnan') std(vel(:,2,2),'omitnan')];
        u_ci = vecnorm(uw_coarse.*uw_std)/u_coarse;
        adcp.burst.u_coarse = u_coarse;
        adcp.burst.u_coarse_ci = u_ci;

        % repeat for attitude data
        % find window data
        idxt_att = adcp_all.attitude.time>=t1 & adcp_all.attitude.time<=t2;
        
        % trim data down to time and range for each window
        adcp.attitude.time = adcp.attitude.time(idxt_att);

        adcp.attitude.pressure = adcp.attitude.pressure(idxt);
        adcp.attitude.heading = adcp.attitude.heading(idxt);
        adcp.attitude.pitch = adcp.attitude.pitch(idxt);
        adcp.attitude.roll = adcp.attitude.roll(idxt);
        adcp.attitude.batteryvoltage = adcp.attitude.batteryvoltage(idxt);
        adcp.attitude.ahrsgyro = adcp.attitude.ahrsgyro(idxt,:);
        adcp.attitude.ahrsrotationmatrix = adcp.attitude.ahrsrotationmatrix(idxt,:);

        adcp.attitude.pitch_ice = adcp.attitude.roll-270;
        adcp.attitude.pitch_ice(adcp.attitude.pitch_ice<-180) = adcp.attitude.pitch_ice(adcp.attitude.pitch_ice<-180) + 360;
        adcp.attitude.roll_ice = -adcp.attitude.pitch;
        adcp.attitude.roll_ice(adcp.attitude.roll_ice<-180) = adcp.attitude.roll_ice(adcp.attitude.roll_ice<-180) + 360;
    
        % save trimmed and processed data structure
        proc_file = fullfile('F:meltstake\data\proc',dep_names{i},sprintf('adcp%d.mat',j));
        fprintf('saving %s\n',proc_file)
        save(proc_file,'adcp')
    
    end
end







