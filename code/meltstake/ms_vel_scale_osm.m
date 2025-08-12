% Script to calculate velocity scales for OSM results. This is a rough
% draft of the scales used for my thesis and melt rate paper.
%
% Instructions:
% set i (deployment) and j (window) values
%
% KJW
% 6 Feb 2023

clear

tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

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

%%
for i = 15%:n_deps

    % time window
    fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    idx_dep = dep_nums(i)==ms_tbl.Number;
    t_start = ms_tbl.Start(idx_dep);
    t_end = ms_tbl.End(idx_dep);
    
    % load data
    load(fullfile('F:meltstake\data\raw',dep_names{i},'adcp','adcp.mat'))
    
    % window
    n_wndws = sum(ms_tbl.Number==dep_nums(i));
    % select window
    for j = 1%:n_wndws
    
        fprintf('--- Window %d ---\n',j)
        
        row_num = find(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        t1 = ms_tbl.Start(row_num);
        t2 = ms_tbl.End(row_num);
        
        % process data
        cor_min = ms_tbl.cor_min(row_num);
    %     fprintf('cor_min: %d\n',cor_min)
        amp_min = 0;

        % unwrap
        adcp.burst.vel(adcp.burst.cor<cor_min) = nan;
        adcp.burst.vel = unwrapVel(adcp.burst.vel,.2,[],32);

        adcp = msADCPTransform(adcp,cor_min,amp_min);
        
        % processing plots
        figure(1001); clf
        adcpQuickPlot(figure(1001),adcp,'vel',0.15*[-1 1],NaT,[0 1],1);
        
        figure(1002); clf
        adcpQuickPlot(figure(1002),adcp,'vel_ice',0.12*[-1 1],NaT,[0 1],1);
        
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
        
        idxr = adcp.burst.range<=r_max(row_num);
        range = adcp.burst.range(idxr);
        
        % mask out data
        vel = adcp.burst.vel_ice(idxt,idxr,1:3);
        vel_allr = adcp.burst.vel_ice(idxt,:,1:3);
        
        % set up second structure for further plotting
        adcp2 = struct('burst',struct);
        adcp2.burst.time = time;
        adcp2.burst.range = range;
        adcp2.burst.vel = vel;
        
        % hanning smooth
        hann_dt = 2; % s
        hann_width = hann_dt*adcp.burst.samplerate;
        for p = 1:size(vel,2)
            for q = 1:size(vel,3)
                if p<=sum(idxr)
                    vel(:,p,q) = hannFilter(vel(:,p,q),hann_width);
                end
                vel_allr(:,p,q) = hannFilter(vel_allr(:,p,q),hann_width);
            end
        end
        adcp2.burst.vel_smth = vel;
        
        % post-mask plots
        figure(1003); clf
        adcpQuickPlot(figure(1003),adcp2,'vel',0.08*[-1 1],NaT,NaN,1);
        
        figure(1004); clf
        adcpQuickPlot(figure(1004),adcp2,'vel_smth',0.08*[-1 1],NaT,NaN,1);
        
        % velocity magnitude
        vel_vol = squeeze(mean(vel,2,'omitnan')); % volume average of u,w,v
        
        vel_mag_prof_max = [max(vecnorm(vel,2,3),[],2) max(vecnorm(vel(:,:,1:2),2,3),[],2)]; % max values of velocity profiles (with v, w/out v)
        
        vel_mag = [vecnorm(vel_vol,2,2) vecnorm(vel_vol(:,1:2),2,2) ... % volume-averaged, then normed (with v, w/out v)
                   mean(vecnorm(vel,2,3),2,'omitnan') mean(vecnorm(vel(:,:,1:2),2,3),2,'omitnan')]; % normed, then volume-averaged
        
        % mean vel profile
        vel_prof = squeeze(mean(vel,1,'omitnan'));
        vel_prof_std = squeeze(std(vel,1,'omitnan'));
        
        vel_prof_rms = sqrt(squeeze(mean(vel.^2,1,'omitnan')));
        vel_prof_rms_max = max(vel_prof_rms,[],1);
        
        % post bin-ave plots
        figure(1005); clf
        subplot(2,1,1)
        plot(time,vel_vol);
        ylabel('vel vol average')
        
        subplot(2,1,2); hold on
        plot(time,vel_mag)
        plot(time,vel_mag_prof_max,'--')
        ylabel('|vel|')
        
        % u scales
        % 1. volume-averaged then normed
        % 2. volume-averaged then normed (no v)
        % 3. normed then volume-averaged
        % 4. normed then volume-averaged (no v)
        % 5. max of magnitude profile
        % 6. max of magnitude profile (no v)
        % 7. max u(t)
        % 8. max w(t)
        % 9. max v(t)
        % 10. max rms(u(y))
        % 11. max rms(w(y))
        % 12. max rms(v(y))
        % 13. max abs(u(y)), retain sign
        % 14. max abs(w(y)), retain sign
        % 15. max abs(v(y)), retain sign
        
        lbls = {'volume-averaged then normed',...
        'volume-averaged then normed (no v)',...
        'normed then volume-averaged',...
        'normed then volume-averaged (no v)',...
        'max of magnitude profile',...
        'max of magnitude profile (no v)',...
        'max(|u(y)|)(t)',...
        'max(|w(y)|)(t)',...
        'max(|v(y)|)(t)',...
        'max rms(u(y))',...
        'max rms(w(y))',...
        'max rms(v(y))'};
        
        u_scales = nan(nt,9);
        u_scales(:,1:4) = vel_mag;
        u_scales(:,5:6) = vel_mag_prof_max;
        u_scales(:,7) = max(abs(vel_allr(:,:,1)),[],2);
        u_scales(:,8) = max(abs(vel_allr(:,:,2)),[],2);
        u_scales(:,9) = max(abs(vel_allr(:,:,3)),[],2);
        
        u_mean = mean(u_scales,'omitnan');
        u_std = std(u_scales,'omitnan');
        
        figure(1016); clf
        for k = 1:9
            subplot(3,3,k); hold on
            plot(time,u_scales)
            plot(time,u_scales(:,k),'k')
            title(lbls{k})
            fprintf('%.3f (%.3f)\n',u_mean(k),u_std(k))
        end
        for k = 1:3
            fprintf('%.3f (0.0)\n',vel_prof_rms_max(k))
        end
        for k = 1:3
            [~,idxy] = max(abs(vel_prof(:,k)));
            fprintf('%.3f (%.3f)\n',vel_prof(idxy,k),vel_prof_std(idxy,k))
        end
        
        figure(1007); clf
        hold on
        plot(range,0*range,'k-')
        for k = 1:3
            errorbar(range,vel_prof(:,k),vel_prof_std(:,k),'color',colors(k),'linewidth',1)
            plot(range,vel_prof(:,k),'ko','markerfacecolor',colors(k))
            plot(range,vel_prof_rms(:,k),'x--','color',colors(k),'linewidth',1)
        end
        title('mean u(y),w(y),v(y)')
        
        % decorrelation scale
        figure(1008); clf
        nk = 400;
        tau = (0:nk)/adcp.burst.samplerate;
        rho = nan(nk+1,3);
        tau_decorr = [nan nan nan];
        for k = 1:3
            rhok = laggedCrossCorr(vel_vol(:,k),vel_vol(:,k),nk);
            rho(:,k) = rhok(nk+1:end);
            if any(rhok<=0.5)
                tau_decorr(k) = tau(find(rho(:,k)<=0.5,1,'first'));
            end
            fprintf('%.1f  ',tau_decorr(k))
        end
        plot(tau,rho,'.-')
        fprintf('[%.1f]\n',mean(tau_decorr,'omitnan'))
    
        % save scales and timestamps
        % save(fullfile('F:meltstake\data\proc',dep_names{i},sprintf('u_scales_%d.mat',j)),'u_scales','time','lbls')
    
    end
end







