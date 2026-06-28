% Script to process open-fjord ctd profiles during meltstake deployments
%
% KJW
% 26 Jun 2026

clear

addpath('..')

proc_dir = 'F:/meltstake/data/proc';

% save output
save_output = 0;

% ms deployment info
ms_tbl = loadMSInfo(26:28,'segments');
[dep_nums,uidx] = unique(ms_tbl.Number);
dep_names = ms_tbl.Folder(uidx);
ndeps = length(dep_nums);
dep_depth = ms_tbl.depth(uidx);

% load ctd data
lgd_lbls = cell(ndeps,1);
ctd = cell(ndeps,1);
for i = 1:ndeps
    load_data = load(fullfile(proc_dir,dep_names{i},'ctd.mat'));
    ctd{i} = load_data.ctd;
    % use polly data if available
    npp = size(ctd{i}(strcmp('polly',{ctd{i}.vessel})).T,2);
    if npp
        ctd{i} = ctd{i}(strcmp('polly',{ctd{i}.vessel}));
    else
        ctd{i} = ctd{i}(~strcmp('polly',{ctd{i}.vessel}));
    end
    t_start = ms_tbl.Start(uidx(i));
    t_start.Format = 'MMMdd yyyy';
    lgd_lbls{i} = sprintf('dep %d (%s)',i,t_start); 
end 
ctd = cell2mat(ctd);

%% data processing
for i = 1:ndeps
    [nz,np] = size(ctd(i).T);
    % clean up ground-strikes in salinity
    iz_min = 0;
    for j = 1:np
        iz_bottom = nz - find(~isnan(flip(ctd(i).S(:,j))),1) + 1;
        [~,idx_filt] = sigmaFilter(ctd(i).S(iz_bottom-10:iz_bottom,j),1,1,2);
        if any(idx_filt(end-1:end))
            % if any outliers detected, clear bottom meter
            ctd(i).S(iz_bottom-3:iz_bottom,j) = ctd(i).S(iz_bottom-4,j);
            ctd(i).PD(iz_bottom-3:iz_bottom,j) = ctd(i).PD(iz_bottom-4,j);
        end
        iz_min = max([iz_min iz_bottom]);
    end

    % calculate mean profiles
    ctd(i).depth_mean = ctd(i).depth(1:iz_min);
    ctd(i).T_mean = mean(ctd(i).T(1:iz_min,:),2,'omitnan');
    ctd(i).S_mean = mean(ctd(i).S(1:iz_min,:),2,'omitnan');
    ctd(i).rho_mean = mean(ctd(i).PD(1:iz_min,:),2,'omitnan')+1000;

    ctd(i).T_std = std(ctd(i).T(1:iz_min,:),0,2,'omitnan');
    ctd(i).T_min = min(ctd(i).T(1:iz_min,:),[],2,'omitnan');
    ctd(i).T_max = max(ctd(i).T(1:iz_min,:),[],2,'omitnan');
    ctd(i).S_std = std(ctd(i).S(1:iz_min,:),0,2,'omitnan');
    ctd(i).S_min = min(ctd(i).S(1:iz_min,:),[],2,'omitnan');
    ctd(i).S_max = max(ctd(i).S(1:iz_min,:),[],2,'omitnan');
    ctd(i).rho_std = std(ctd(i).PD(1:iz_min,:),0,2,'omitnan');
    ctd(i).rho_min = min(ctd(i).PD(1:iz_min,:),[],2,'omitnan')+1000;
    ctd(i).rho_max = max(ctd(i).PD(1:iz_min,:),[],2,'omitnan')+1000;

    % calculate N2
    dz = diff(ctd(i).depth(1:2));
    rho_mean_smth = hannFilter(ctd(i).rho_mean,round(3/dz)+1);
    drho_dz = gradient(rho_mean_smth)/dz;
    ctd(i).N2 = 9.81*drho_dz./rho_mean_smth/(2*pi)^2;
    ctd(i).N = ctd(i).N2; % lol
    ctd(i).N(ctd(i).N<0) = nan;
    ctd(i).N = sqrt(ctd(i).N);

    % find melstake depth index
    [~,ctd(i).zidx] = min(abs(ctd(i).depth_mean-dep_depth(i)));

    % find range of N values
    z_range = 5; % m
    dzi = round(z_range/dz);
    ctd(i).N_ms = extrema(ctd(i).N((ctd(i).zidx-dzi):(ctd(i).zidx+dzi)));

end

%% save
if save_output
    save(fullfile(proc_dir,'glacier_2024_ctd.mat'),'ctd')
end