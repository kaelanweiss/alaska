% Script to calculate velocity scales for meltstake deployments. This is a
% more finished version of the one created for OSM. The choice of maximum
% of the magnitude profile mimics a freestream velocity.
%
% Instructions:
% set i (deployment) and j (window) values
%
% KJW
% 4 Apr 2024
%
% NOTE: for quick calculation, I added/changed lines 22, 42, 66

clear

ms_tbl = loadMSInfo(26:28,'manualwindows');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');
r_max = ms_tbl.rmax;
n_deps = length(dep_nums);

% times to exclude due to bad data, corresponding deployment number
t_exclude = [datetime(2023,5,29,21,59,0) datetime(2023,5,29,22,03,0);...
             datetime(2023,5,30,1,41,0) datetime(2023,5,30,1,47,0);...
             datetime(2023,7,9,21,30,0) datetime(2023,7,9,21,40,0);...
             datetime(2023,7,9,21,45,0) datetime(2023,7,9,21,51,30);...
             datetime(2023,7,9,23,22,0) datetime(2023,7,9,23,28,0)];
t_exclude_dep = [2;...
                 3;...
                 8;...
                 8;...
                 8];

u_max = nan(32,2);
up2 = nan(32,1);
v_rms = nan(32,1);
u_coarse = nan(32,2);
pos = 1;

%%
for i = 1:n_deps

    % time window
    fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    idx_dep = dep_nums(i)==ms_tbl.Number;
    t_start = ms_tbl.Start(idx_dep);
    t_end = ms_tbl.End(idx_dep);
    
    % load data (this should be switched to processed ADCP files)
    load(fullfile('F:meltstake\data\raw',dep_names{i},'adcp','adcp.mat'))
    adcp_load = adcp;
    
    % window
    n_wndws = sum(ms_tbl.Number==dep_nums(i));
    % select window
    for j = 1:n_wndws
    
        fprintf('--- Window %d ---\n',j)
        
        row_num = find(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        t1 = ms_tbl.Start(row_num);
        t2 = ms_tbl.End(row_num);
        
        % processing values
        cor_min = ms_tbl.cor_min(row_num);
        amp_min = 0;

        % unwrap
        adcp = adcp_load;
        adcp.burst.vel(adcp.burst.cor<cor_min) = nan;
%         adcp.burst.vel = unwrapVel(adcp.burst.vel,.2,[],32);

        % smooth across 3 points in the profile
        for k = 1:size(adcp.burst.vel,1)
            for m = 1:size(adcp.burst.vel,3)
                adcp.burst.vel(k,:,m) = hannFilter(squeeze(adcp.burst.vel(k,:,m)),3);
            end
        end
        
        % transform velocity
        adcp = msADCPTransform(adcp,cor_min,amp_min);
        
        % processing plots
%         figure(1001); clf
%         adcpQuickPlot(figure(1001),adcp,'vel',extrema(adcp.burst.vel),NaT,[0 1],1);
%         
%         figure(1002); clf
%         adcpQuickPlot(figure(1002),adcp,'vel_ice',0.12*[-1 1],NaT,[0 1],1);
        
        % find window data
        time = dn2dt(adcp.burst.time);
        idxt = time>=t1 & time<=t2;
        % exclude specified times
        for k = 1:size(t_exclude,1)
            if t_exclude_dep(k)==dep_nums(i)
                idxk = time>=t_exclude(k,1) & time<=t_exclude(k,2);
                idxt = idxt & ~idxk;
            end
        end
        time = time(idxt);
        nt = length(time);

        % range of good data
        idxr = adcp.burst.range<=r_max(row_num);
        range = adcp.burst.range(idxr);
        
        % mask out data
        vel = adcp.burst.vel_ice(idxt,idxr,1:3);
        
        % set up second structure for further plotting
        adcp2 = struct('burst',struct);
        adcp2.burst.time = time;
        adcp2.burst.range = range;
        adcp2.burst.vel = vel;
        
        % hanning smooth
        hann_dt = 1; % s
        hann_width = hann_dt*adcp.burst.samplerate;
        for p = 1:size(vel,2)
            for q = 1:size(vel,3)
%                 if p<=sum(idxr)
                    vel(:,p,q) = hannFilter(vel(:,p,q),hann_width);
%                 end
            end
        end
        adcp2.burst.vel_smth = vel;
        
        % post-mask plots
%         figure(1003); clf
%         adcpQuickPlot(figure(1003),adcp2,'vel',0.08*[-1 1],NaT,NaN,1);
%         
%         figure(1004); clf
%         adcpQuickPlot(figure(1004),adcp2,'vel_smth',0.08*[-1 1],NaT,NaN,1);
        
        % u'^2 (energy)
        % mean profiles of u,w,v
        vel_prof = mean(vel(:,:,1:3),1,'omitnan');
        % u'
        vel_prime = vel - repmat(vel_prof,[size(vel,1) 1 1]);
        vel_prime2 = vel_prime.^2;
        vel_prime2_mean = squeeze(mean(mean(vel_prime2,1,'omitnan'),2,'omitnan'));
        uprime2 = sum(vel_prime2_mean);
        
        % magnitude of u,w
        vel_mag = vecnorm(vel(:,:,1:2),2,3);

        % profile max
        vel_max = max(vel_mag,[],2);
        vel_max_smth = hannFilter(vel_max,4*5*60);
        vel_max_prime = vel_max - vel_max_smth;
        vel_rms = sqrt(mean(vel_max_prime.^2,'omitnan'));

        % "outer" scale (time mean of u, w next to ADCP)
        uw_outer = mean(squeeze(vel(:,2,1:2)),'omitnan');
        U_outer = vecnorm(uw_outer);
        uw_std = [std(vel(:,2,1),'omitnan') std(vel(:,2,2),'omitnan')];
        U_ci = vecnorm(uw_outer.*uw_std)/U_outer;

%         figure(pos); clf; hold on
%         plot(vel_max,'k')
%         plot(vel_max_smth)
%         plot(vel_max_prime)
%         plot(vel_rms+0*vel_max)
        
        u_max(pos,:) = [mean(vel_max,'omitnan') std(vel_max,'omitnan')];
        up2(pos) = uprime2;
        v_rms(pos) = vel_rms;
        u_coarse(pos,:) = [U_outer U_ci];

        fprintf('%d %.3f (%.3f) %.6f %.3f %.3f\n',pos,mean(vel_max,'omitnan'),std(vel_max,'omitnan'),uprime2,vel_rms,U_outer)
        
        pos = pos+1;
    
    end
end







