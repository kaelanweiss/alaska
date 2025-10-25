% Script to perform bin mapping of the ADCP data relative to the ice
% locations (determined using script "manual_adcp_ice_edges.m").
%
% KJW
% 23 Oct 2025

clear

raw_dir = 'F:meltstake/data/raw';
proc_dir = 'F:meltstake/data/proc';

% choose deployment
dep_num = 28;
seg_tbl = loadMSInfo(dep_num,'segments');
ns = size(seg_tbl,1);

% load data
load(fullfile(raw_dir,seg_tbl.Folder{1},'adcp','adcp.mat'))
load(fullfile(proc_dir,seg_tbl.Folder{1},'adcp_ice_edges.mat'))
adcp.burst.range = cast(adcp.burst.range,'double');
%% plot raw data with edges for reference
fig1 = figure(1); clf
ax1 = adcpQuickPlot(fig1,adcp,'vel',0.1*[-1 1],NaT,[0 max(seg_tbl.rmax)+0.3],1);
for i = 1:length(ax1)
    hold(ax1(i),'on')
end

for i = 1:ns
    idxt = adcp.burst.time>=seg_tbl.Start(i) & adcp.burst.time<=seg_tbl.End(i);
    ti = adcp.burst.time(idxt);

    for j = 1:adcp.burst.nbeams
        tg = datetime(adcp_ice_edges{i}{j}(:,1),'convertfrom','datenum'); % times
        rg = adcp_ice_edges{i}{j}(:,2); % ranges
        if length(tg) == 1
            tg = tg+[0; seconds(0.1)];
            rg = rg*[1;1];
        end
        icepos = interp1(tg,rg,ti,'nearest','extrap'); % interpolated ice position
        plot(ax1(j),ti,icepos+adcp.burst.cellsize/2,'g-')
    end
end

%% remap and transform
% remapping distance
z0 = nan(ns,1);
z1 = adcp.burst.range(1);

for i = 1:ns
    % time axis
    idxt = adcp.burst.time>=seg_tbl.Start(i) & adcp.burst.time<=seg_tbl.End(i);
    ti = adcp.burst.time(idxt);
    nt = length(ti);

    % build ice edge time series out of manual points
    ice_edge = nan(length(ti),adcp.burst.nbeams);
    for j = 1:adcp.burst.nbeams
        tg = datetime(adcp_ice_edges{i}{j}(:,1),'convertfrom','datenum'); % times
        rg = adcp_ice_edges{i}{j}(:,2); % ranges
        if length(tg) == 1 % duplicate if only one point exists
            tg = tg+[0; seconds(0.125)];
            rg = rg*[1;1];
        end
        ice_edge(:,j) = interp1(tg,rg,ti,'nearest','extrap'); % interpolated ice position
    end

    % remapping distance is the mean ice position
    z0(i) = mean(ice_edge,'all');

    % remapping axis
    nz = ceil(2*(z0(i)-z1)/adcp.burst.cellsize); % aim for twice the original spatial resolution
    zq = linspace(z1,z0(i),nz);
    vel_rmp = nan(nt,nz,adcp.burst.nbeams);

    for j = 1:adcp.burst.nbeams        
        % remap velocities (keep first adcp bin the same, interpolate
        if isfield(adcp.burst,'vel_unw')
            velj = adcp.burst.vel_unw(idxt,:,j);
        else
            velj = adcp.burst.vel(idxt,:,j);
        end
        % QC
        velj(adcp.burst.cor(idxt,:,j)<adcp.burst.processing.cor_min) = nan;

        for k = 1:nt
            % compute scaled bin positions
            zw = ice_edge(k,j);
            zk = (z0(i)-z1)/(zw-z1)*(adcp.burst.range-z1) + z1;
            % interpolate
            velk = interp1(zk,velj(k,:),zq,'linear');
            vel_rmp(k,:,j) = velk;
        end

    end

    % angles for rotating velocity, right handed along instrument x (ms right) and y (down) axes, respectively
    angle_horz = adcpIceSlope(ice_edge(:,3),ice_edge(:,1));
    angle_vert = adcpIceSlope(ice_edge(:,4),ice_edge(:,2));

    % load processed adcp
    adcp_rmp = load(fullfile(proc_dir,seg_tbl.Folder{1},sprintf('adcp%d.mat',i)));
    adcp_rmp = adcp_rmp.adcp;
    if isfield(adcp_rmp.burst,'vel_unw')
        adcp_rmp.burst = rmfield(adcp_rmp.burst,{'vel_unw'});
    end
    % insert remapped data and extra info
    adcp_rmp.burst.vel = vel_rmp;
    adcp_rmp.burst.cor = 0*vel_rmp+100; % fake, already QCed above
    adcp_rmp.burst.amp = 0*vel_rmp+max(adcp_rmp.burst.amp,[],'all'); % also fake
    adcp_rmp.burst.range = zq;
    adcp_rmp.burst.ice_edge = ice_edge;
    adcp_rmp.burst.angle_x_rot = angle_vert;
    adcp_rmp.burst.angle_y_rot = angle_horz;

    % perform regular transform
    adcp_rmp = msADCPTransform(adcp_rmp);
    adcp_rmp.burst.vel_flat = adcp_rmp.burst.vel_ice;
    adcp_rmp.burst.vel_ice = nan*adcp_rmp.burst.vel_ice;

    % perform velocity rotations according to ice distances
    phi_x = mean(angle_vert);
    phi_y = mean(angle_horz);
    cx = cos(phi_x); sx = sin(phi_x);
    cy = cos(phi_y); sy = sin(phi_y);
    % rotation matrix
    %      (u   w  v5 v1-4  err)
    Rx = [  1   0   0   0   0;... % u
            0  cx  sx   0   0;... % w'
            0 -sx  cx   0   0;... % v5'
            0 -sx   0  cx   0;... % v1-4'
            0   0   0   0   1];   % err
    
    %      (u   w  v5 v1-4  err)
    Ry = [  cy   0 -sy   0   0;... % u'
             0   1   0   0   0;... % w'
            sy   0  cy   0   0;... % v5
            sy   0   0  cy   0;... % v1-4
             0   0   0   0   1];   % err
    R = Rx*Ry;
    % loop through cells
    for j = 1:nz
        velj = squeeze(adcp_rmp.burst.vel_flat(:,j,:))';
        velj = R*velj;
        adcp_rmp.burst.vel_ice(:,j,:) = permute(velj,[2 3 1]);
    end

    % save output
    adcp_raw = adcp;
    adcp = adcp_rmp;
    save(fullfile(proc_dir,seg_tbl.Folder{1},'adcp_bin_map',sprintf('adcp%d.mat',i)),'adcp')
    adcp = adcp_raw;
end


