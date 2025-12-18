% Script to calculate PSDs from ADCP observations at the glacier. These
% data will be used in other figures, mainly figure 3.
%
% KJW
% 12 Nov 2025
clear

raw_dir = 'F:/meltstake/data/raw';

% collect velocity data
% meltstake
fprintf('Meltstake data...\n')
dep_tbl = loadMSInfo(26:28);
seg_tbl = loadMSInfo(26:28,'segments');
dep_nums = dep_tbl.Number;
dep_names = dep_tbl.Folder;
ndeps = length(dep_nums);

adcp(ndeps) = struct('dep_name',[],'fs',[],'time',[],'range',[],'vel',[],'vel_ice',[]);
for i = 1:ndeps
    % loading
    adcp_load = load(fullfile(raw_dir,dep_names{i},'adcp','adcp.mat'));
    adcp_i = adcp_load.adcp;
    % indexing
    idxd = strcmp(seg_tbl.Folder,dep_names{i});
    t1 = min(seg_tbl.Start(idxd));
    t2 = max(seg_tbl.End(idxd));
    idxt = adcp_i.burst.time >= t1 & adcp_i.burst.time <= t2;
    idxr = adcp_i.burst.range <= min(seg_tbl.rmax(idxd));
    % grabbing data
    adcp(i).dep_name = dep_names{i};
    adcp(i).fs = adcp_i.burst.samplerate;
    adcp(i).time = adcp_i.burst.time(idxt);
    adcp(i).range = adcp_i.burst.range(idxr);
    adcp(i).vel = adcp_i.burst.vel(idxt,idxr,:);
    adcp(i).vel_ice = adcp_i.burst.vel_ice(idxt,idxr,[1 3 2]);
end
clear adcp_load adcp_i

%% calculate psds 
for i = 1:ndeps
    % preallocate psds
    [nt,nc,nb] = size(adcp(i).vel);
    nf = ceil((nt+2)/2);
    
    P_i = nan(nf,nc,nb);
    P_ice_i = nan(nf,nc,3);
    
    % beam velocities
    v_mean = mean(adcp(i).vel,1,'omitnan');
    for j = 1:nb
        vj = adcp(i).vel(:,:,j)-repmat(v_mean(:,:,j),[nt 1]);
        P_i(:,:,j) = psd_matlab(vj,adcp(i).fs,'notaper');
    end

    % ice velocities
    % interpolate over nans from the processing routine
    for j = 1:nc
        idx_nan = isnan(adcp(i).vel_ice(:,j,1));
        for k = 1:3
            adcp(i).vel_ice(idx_nan,j,k) = interp1(adcp(i).time(~idx_nan),adcp(i).vel_ice(~idx_nan,j,k),adcp(i).time(idx_nan),'linear','extrap');
        end
    end
    % psd
    v_ice_mean = mean(adcp(i).vel_ice,1,'omitnan');
    for j = 1:3
        vj = adcp(i).vel_ice(:,:,j)-repmat(v_ice_mean(:,:,j),[nt 1]);
        [P_ice_i(:,:,j),fi] = psd_matlab(vj,adcp(i).fs,'notaper');
    end

    % store values
    adcp(i).f = fi;
    adcp(i).P = P_i;
    adcp(i).P_ice = P_ice_i;
end

%% save
% save platform_comparison_data\ms_psd.mat adcp




