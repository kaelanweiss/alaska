% Script to calculate melt rates from observations of u, T, and S according
% to the 3-eqn model, K&M15, Schulz22, and FitzMaurice17.
%
% KJW
% 7 Apr 2024

clear

ms_tbl = loadMSInfo(26:28,'manualwindows');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');
n_deps = length(dep_nums);

% data
proc_path = 'F:/meltstake/data/proc';

% set smoothing width (Hanning)
hann_int = .5; % s

m_array = cell(100,1);
pos = 1;

% Schulz threshold velocity
u_thresh = 0.05;

m_desc = {'3equation','KerrMcConnochie2015convective',sprintf('3eqnSchulz22(%.3fm/s)',u_thresh),'FitzMaurice2017bulk',...
    '3eqn_coarse',sprintf('3eqnSchulz22_coarse(%.3fm/s)',u_thresh),'3eqn_noS'};

% loop through deployments
for i = 1:n_deps

    fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    idx_dep = dep_nums(i)==ms_tbl.Number;
    
    % loop through sections
    n_wndws = sum(ms_tbl.Number==dep_nums(i));
    for j = 1:n_wndws

        fprintf('window %d\n',j)

        % load data
        load(fullfile(proc_path,dep_names{i},['adcp' num2str(j) '.mat']))
        load(fullfile(proc_path,dep_names{i},['T' num2str(j) '.mat']))
%         load(fullfile(proc_path,dep_names{i},['ctd' num2str(j) '.mat']))

        % calculate constant salinity scale
%         S = mean(ctd.S,'all','omitnan');
        S = msTable2Vector(ms_tbl.S(idx_dep & ms_tbl.Window == j));

        % smooth U
        vel = adcp.burst.vel_ice(:,:,1:2); % u,w
        % smooth across 3 points in the profile
        for k = 1:size(vel,1)
            for m = 1:size(vel,3)
                vel(k,:,m) = hannFilter(squeeze(vel(k,:,m)),3);
            end
        end
        % calculate magnitude
        vel_mag = vecnorm(vel,2,3);
        % profile max and mean
        vel_max = max(vel_mag,[],2);
        u_max = max(vel(:,:,1),[],2);
        w_max = max(vel(:,:,2),[],2);

        % extract subsampled time series according to smoothing window
        t_scale = T.time;
        T_scale = mean(T.T,2,'omitnan');
        U_scale = interp1(adcp.burst.time,vel_max,t_scale,'nearest',nan);
        u_scale = interp1(adcp.burst.time,u_max,t_scale,'nearest',nan);
        w_scale = interp1(adcp.burst.time,w_max,t_scale,'nearest',nan);
        U_coarse = msTable2Vector(ms_tbl.u_coarse(pos));
        T_coarse = msTable2Vector(ms_tbl.T(pos));

        % compute melt rate on the subsampled time axis
        S_scale = S*ones(size(T_scale));

        % 3-equation
        [m1,Tb,Sb] = solve3Eqn(U_scale,T_scale,S_scale,0);
        m1 = m1*86400; % m/s --> m/day

        % KM15
        m2 = (2.5e-7)*(T_scale-Tb).^(4/3)*86400;
%         m2 = (2.5e-7)*(T_scale).^(4/3)*86400;

        % Schulz22
        m3 = solve3EqnSchulz22(U_scale,T_scale,S_scale,0,u_thresh);
        m3 = m3*86400; % m/s --> m/day

        % FitzMaurice17
        L = 10; % m
        Ti = 0; % C
        m4 = solveFM17(U_scale,u_scale,w_scale,T_scale,Ti,L);

        % 3-equation with coarse scales
        m5 = solve3Eqn(U_coarse,T_coarse,S_scale(1),0);
        m5 = m5*86400; % m/s --> m/day

        % Schulz22 with coarse scales
        m6 = solve3EqnSchulz22(U_coarse,T_coarse,S_scale(1),0,u_thresh);
        m6 = m6*86400; % m/s --> m/day

        % 3-equation with S = 0
        [m7,Tb7,Sb7] = solve3Eqn(U_scale,T_scale,0*S_scale,0);
        m7 = m7*86400; % m/s --> m/day
        
        % "save" output
        m_array{pos} = [m1 m2 m3 m4 m5*ones(size(m1)) m6*ones(size(m1)) m7];
        pos = pos + 1;
    end
end
m_array(pos:end) = []; % trim unused space

% calculate mean melt for each
m_mean = nan(length(m_array),size(m_array{1},2));
m_std = nan(length(m_array),size(m_array{1},2));
for i = 1:length(m_array)
    m_mean(i,:) = mean(m_array{i},'omitnan');
    m_std(i,:) = std(m_array{i},'omitnan');
end

%%
save('F:meltstake/data/proc/melt_models_glacier.mat','m_array','m_mean','m_std','m_desc')
save('../../data/melt_models_glacier.mat','m_array','m_mean','m_std','m_desc')

% melt_paper_figure_3
