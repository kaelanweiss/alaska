% Script to calculate melt rate and scales for OSM scatter plot. This
% relies on the manually specified windows in the
% meltstake_deployments.xlsx spreadsheet.
%
% KJW
% 5 Feb 2024

clear

t_pad = minutes(15);

tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');
n_deps = length(dep_names);

% loop through deployments
for i = 1:n_deps
%     fprintf('====== DEPLOYMENT %d (%s) ======\n',dep_nums(i),dep_names{i})
    idx_dep = dep_nums(i)==ms_tbl.Number;
    t_start = ms_tbl.Start(idx_dep);
    t_end = ms_tbl.End(idx_dep);

    % load data
    raw_fldr = fullfile('F:meltstake\data\raw',dep_names{i});
    proc_fldr = fullfile('F:meltstake\data\proc',dep_names{i});
    % Temp
    if exist(fullfile(raw_fldr,'rbr','T.mat'),'file')
        load(fullfile(raw_fldr,'rbr','T.mat'))
        includeT = 1;
        T_all = repmat(T(1).values,[1 length(T)]);
        for j = 2:length(T)
            T_all(:,j) = T(j).values;
        end
    else
        includeT = 0;
    end

    % ADCP
    load(fullfile(raw_fldr,'adcp','adcp.mat'))
    adcp_depth = adcp.attitude.pressure-1;
    cor_min = 25;
    adcp = msADCPTransform(adcp);

    % Sal
    if exist(fullfile(proc_fldr,'ctd.mat'),'file')
        load(fullfile(proc_fldr,'ctd.mat'))
        includeS = 1;
    end

    % loop through windows in the deployment
    nj = sum(idx_dep);
    for j = 1:nj
%         fprintf('--- Window %d ---\n',j)
        t1 = t_start(j);
        t2 = t_end(j);
        
        % window depth
        zj = mean(adcp_depth(dn2dt(adcp.attitude.time)>=t1 & dn2dt(adcp.attitude.time)<=t2));
%         fprintf('Depth=%.1f',zj)
        
        % Temperature
        if includeT
            idxT = dn2dt(T(1).time)>=t1 & dn2dt(T(1).time)<=t2;
            Tj = T_all(idxT,:);
            T_scale = mean(Tj,'all','omitnan');
%             T_scale2 = max(mean(Tj,'omitnan'));
            T_std = std(Tj(:),'omitnan');
%             fprintf('\tT=%.1f (%.1f)',T_scale,T_std)
            
%             figure; hold on
%             plot(dn2dt(T(1).time(idxT)),Tj)
%             errorbar(mean(dn2dt(T(1).time(idxT))),T_scale,T_std,'k','linewidth',1)
%             plot(mean(dn2dt(T(1).time(idxT))),T_scale,'ko','markerfacecolor','k')
%             title(sprintf('Dep %d window %d',dep_nums(i),j))
        end

        % Salinity and CTD T
        if includeS
            idxS = ctd.t0>=(t1-t_pad) & ctd.t0<=(t2+t_pad);
            idxz = abs(ctd.z-round(2*zj)/2)<=1;
            if sum(idxS)==0
                % find  preceding/following profiles
                S_scale = -9;
                S_std = -9;
                Sj = nan*ctd.z(idxz);

                Tctd_scale = -9;
                Tctd_std = -9;
                Tcdtj = nan*ctd.z(idxz);
            else
                Sj = ctd.S(idxz,idxS);
                S_scale = mean(Sj,'all','omitnan');
                S_std = std(Sj(:),'omitnan');
                fprintf('\tS=%.1f (%.1f)',S_scale,S_std)

                Tctdj = ctd.T(idxz,idxS);
                Tctd_scale = mean(Tctdj,'all','omitnan');
                Tctd_std = std(Tctdj(:),'omitnan');
%                 fprintf('\tT_ctd=%.1f (%.1f)',Tctd_scale,Tctd_std)
%                 fprintf('%.1f (%.1f)',Tctd_scale,Tctd_std)
            end
            
%             figure; hold on
%             plot(ctd.z(idxz),Sj)
%             errorbar(zj,S_scale,S_std,'k','linewidth',1)
%             plot(zj,S_scale,'ko','markerfacecolor','k')
%             title(sprintf('Dep %d window %d',dep_nums(i),j))
        else
            fprintf('-9 (-9)\n')
        end

        % ADCP
%         figure(1001); clf
%         adcpQuickPlot(figure(1001),adcp,'vel',0.05*[-1 1],[t1 t2],[0 1],1);
% 
%         figure(1002); clf
%         adcpQuickPlot(figure(1002),adcp,'vel_ice',0.08*[-1 1],[t1 t2],[0 1],1);


        fprintf('\n')

%         input('')
    end
%     fprintf('\n')
end

