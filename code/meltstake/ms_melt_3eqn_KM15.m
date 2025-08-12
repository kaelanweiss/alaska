% Script to calculate melt rates from observations of u, T, and S according
% to the 3-eqn model and K&M15

clear

% tbl_path = 'F:meltstake\metadata\meltstake_deployments.xlsx';
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');
n_deps = length(dep_nums);

% data
proc_path = 'F:/meltstake/data/proc';

% set smoothing width (Hanning)
hann_int = .5; % s

m_all = nan(100,5);

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
        load(fullfile(proc_path,dep_names{i},['ctd' num2str(j) '.mat']))

        % calculate constant salinity scale
        S = mean(ctd.S,'all','omitnan');
        S_std = std(ctd.S(:),'omitnan');

        % smooth T and U
        k_T = round(hann_int/seconds(diff(T.time(1:2))));
        k_U = hann_int*adcp.burst.samplerate;
        % T
%         T_smth = hannFilter(mean(T.T,2,'omitnan'),k_T);
        T_smth = mean(T.T,2,'omitnan');
        % U
        vel = adcp.burst.vel_ice(:,:,1:3); % u,w,v
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
        vel_mean = mean(vel_mag,2,'omitnan');
        % smooth
%         vel_max_smth = hannFilter(vel_max,k_U);
%         vel_mean_smth = hannFilter(vel_mean,k_U);
        vel_max_smth = vel_max;
        vel_mean_smth = vel_mean;

        % extract subsampled time series according to smoothing window
        t_scale = T.time(ceil(k_T/2):k_T:end);
        T_scale = T_smth(ceil(k_T/2):k_T:end);
        U1_scale = interp1(adcp.burst.time,vel_max_smth,t_scale,'nearest',nan);
        U2_scale = interp1(adcp.burst.time,vel_mean_smth,t_scale,'nearest',nan);


        % compute melt rate on the subsampled time axis
        S_scale = S*ones(size(T_scale));
        % 3-equation
        [m1,Tb,Sb] = solve3Eqn(U1_scale,T_scale,S_scale,0);
        m1 = m1*86400; % m/s --> m/day
        % KM15
        m2 = (2.5e-7)*(T_scale-Tb).^(4/3)*86400;

        m3 = max(m1,m2);

        % plot time series of U, T, and m; report mean m with uncertainty
        figure(10*dep_nums(i)+j); clf
        clear ax
        % T
        ax(1) = subplot(3,1,1);
        plot(t_scale,T_scale,'.-','color',colors(2))
        % U
        ax(2) = subplot(3,1,2);
        hold on
        plot(t_scale,U1_scale,'.-','color',colors(1))
        plot(t_scale,U2_scale,'.-','color',colors(6))
        % m
        ax(3) = subplot(3,1,3);
        hold on
        plot(t_scale,m1,'.-','color','k')
        plot(t_scale,m2,'.-','color',0.3*[1 1 1])
        fprintf('%.2f (%.2f) m/day\n%.2f (%.2f) m/day\n%.2f m/day\n',mean(m1,'omitnan'),std(m1,'omitnan'),mean(m2,'omitnan'),std(m2,'omitnan'),mean(m3,'omitnan'))
        
        % "save" output
        m_all(find(isnan(m_all(:,1)),1),:) = [mean(m1,'omitnan'),std(m1,'omitnan'),mean(m2,'omitnan'),std(m2,'omitnan'),mean(m3,'omitnan')];
    end

end

m_all = m_all(~isnan(m_all(:,1)),:);

m_obs = msTable2Vector(ms_tbl.m);
m_obs = m_obs*.24;

figure
for i = 1:3
    subplot(1,3,i)
    hold on; box on
    plot([0 1.5],[0 1.5],'k-')
    plot([0 1.5],3*[0 1.5],'k-')
    plot(m_all(:,2*i-1),m_obs,'ko','markersize',4,'markerfacecolor',colors(1))
    xlim([0 1.5])
    ylim([0 2.5])
    if i == 1
        title(sprintf('%.2fs',hann_int))
    end
end
