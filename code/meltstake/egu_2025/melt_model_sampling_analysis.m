% script to run sample rate analysis on glacier data for EGU 2025
%
% KJW
% 14 Apr 2025

clear

addpath('..')

ms_tbl = loadMSInfo(26:28,'manualwindows');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');
n_deps = length(dep_nums);

% data
proc_path = 'F:/meltstake/data/proc';

% set smoothing width (Hanning)
tau = 2.^(-1:10); % s

% preallocate
m_mean = nan(size(ms_tbl,1),length(tau)+1,2);

m_desc = {'3equation','3eqn_coarse'};

% depug plotting
fig = figure(99); clf
ax = axes(fig);

% loop through deployments
for i = 1:n_deps

    fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    idx_dep = dep_nums(i)==ms_tbl.Number;
    
    % loop through sections
    n_wndws = sum(ms_tbl.Number==dep_nums(i));
    for j = 1:n_wndws

        fprintf('window %d\n',j)

        % row number
        row_num = ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j;

        % load data
        load(fullfile(proc_path,dep_names{i},['adcp' num2str(j) '.mat']))
        load(fullfile(proc_path,dep_names{i},['T' num2str(j) '.mat']))

        % calculate constant salinity scale
        S = msTable2Vector(ms_tbl.S(idx_dep & ms_tbl.Window == j));

        % velocity
        % profile max and mean
        u_max = adcp.burst.u_max;
        % component timeseries
%         u = mean(adcp.burst.vel_ice(:,:,1),2,'omitnan');
%         w = mean(adcp.burst.vel_ice(:,:,2),2,'omitnan');
        u = adcp.burst.vel_ice(:,2,1);
        w = adcp.burst.vel_ice(:,2,2);

        % extract subsampled time series according to smoothing window
        t_scale = T.time;
        T_scale = mean(T.T,2,'omitnan');
        U_max = interp1(adcp.burst.time,u_max,t_scale,'nearest',nan);
        U = interp1(adcp.burst.time,u,t_scale,'nearest',nan);
        W = interp1(adcp.burst.time,u,t_scale,'nearest',nan);
        vel_ice = nan([length(t_scale) size(adcp.burst.vel,2) 2]);
        for p = 1:size(vel_ice,2)
            for q = 1:2
                vel_ice(:,p,q) = interp1(adcp.burst.time,adcp.burst.vel_ice(:,p,q),t_scale,'nearest',nan);
            end
        end

        U_coarse = adcp.burst.u_coarse;
        T_coarse = msTable2Vector(ms_tbl.T(row_num));

        % compute melt rate on the subsampled time axis
        S_scale = S*ones(size(T_scale));

        % debug plotting
        cla(ax)
        hold(ax,'on')
        plot(t_scale,vecnorm([U W],2,2),'k-')
        plot(extrema(t_scale),U_coarse*[1 1],'r-')
        for k = 1:length(tau)
            % smooth inputs
%             Umaxk = hannFilter(U_max,2*tau(k)); % u_max
            Umaxk = calculateUMax(vel_ice,2*tau(k));
            Tk = hannFilter(T_scale,2*tau(k));
            uk = hannFilter(U,2*tau(k));
            wk = hannFilter(W,2*tau(k));
            Umagk = vecnorm([uk wk],2,2);
            % solve 3eqn
            m_max = solve3Eqn(Umaxk,Tk,S_scale,0);
            m_crs = solve3Eqn(Umagk,Tk,S_scale,0);
            % store mean melt rates
            m_mean(row_num,k,1) = mean(m_max,'omitnan')*86400;
            m_mean(row_num,k,2) = mean(m_crs,'omitnan')*86400;

            % debug plot
            if k < length(tau)
                plot(t_scale,Umagk,'k-')
            else
                plot(t_scale,Umagk,'r-')
            end
            
        end
        m_mean(row_num,end,1) = solve3Eqn(mean(U_max,'omitnan'),T_coarse,S_scale(1),0)*86400;
        m_mean(row_num,end,2) = solve3Eqn(U_coarse,T_coarse,S_scale(1),0)*86400;
    end
end


%%
% save melt_3eqn_tau m_mean tau
% save('F:meltstake/data/proc/melt_models_glacier.mat','m_array','m_mean','m_std','m_desc')
% save('../../data/melt_models_glacier.mat','m_array','m_mean','m_std','m_desc')

% melt_paper_figure_3

function u_max = calculateUMax(vel,hann_width)
% function to calculate u_max from the 3D vel input
    [nt,nc,nb] = size(vel);
    
    % smooth components in time
    for i = 1:nc
        for j = 1:nb
            vel(:,i,j) = hannFilter(vel(:,i,j),hann_width);
        end
    end
    
    % smooth along profiles
    for i = 1:nt
        for j = 1:nb
            vel(i,:,j) = hannFilter(vel(i,:,j),3);
        end
    end

    % velocity magnitude
    vel_mag = vecnorm(vel,2,3);
    % u_max
    u_max = max(vel_mag,[],2);

end