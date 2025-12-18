% Script to plot plan view distributions of ocean velocity from Polly/Aries 
% to contextualize MS observations.
%
% KJW
% 10 Dec 2025
clear

gd_proc_dir = 'G:/Shared drives/Ice-ocean-interactions/fieldwork_docs_and_data/LeConte2406/data/processed';

polly_dep = {'20240712_162759','20240715_200834','20240716_211642'};
aries_dep = {'20240712_164057','20240715_163629','20240716_191440'};
term_dep = {{'07122024','01'},{'07152024','01'},{'07162024','01'}};
ms_depth = [21 43 52];

R_e = 6.4e6; % [m] radius of earth
%% load
grid_vel(3) = struct('term',[],'ms_loc',[],'lat',[],'lon',[],'vel',[],'n',[],'depth',[]);
for i = 1:3
% i = 3;

    polly = getfield(load(fullfile(gd_proc_dir,'Polly',sprintf('deploy_%s',polly_dep{i}),sprintf('UBOX02_deploy_%s',polly_dep{i}),sprintf('adcp_deploy_%s.mat',polly_dep{i}))),'adcp');
    aries = getfield(load(fullfile(gd_proc_dir,'Aries',sprintf('deploy_%s',aries_dep{i}),sprintf('UBOX04_deploy_%s',aries_dep{i}),sprintf('adcp_deploy_%s.mat',aries_dep{i}))),'adcp');
    term = load(fullfile(gd_proc_dir,'Drone/LeconteTerminus',term_dep{i}{1},term_dep{i}{2},'TermLatLon',sprintf('LeconteTermVert_%s_%d.mat',term_dep{i}{1},str2double(term_dep{i}{2}))));
    % ctd = getfield(load(fullfile(gd_proc_dir,'Polly',sprintf('deploy_%s',polly_dep{i}),'CTD',sprintf('adcp_deploy_%s.mat',polly_dep{i}))),'adcp');
    
    dep_tbl = loadMSInfo(25+i);
    ms_loc = str2double(strsplit(dep_tbl.Location{1},', '));
    
    polly.nuc_time = datetime(polly.nuc_time,'convertfrom','datenum');
    aries.nuc_time = datetime(aries.nuc_time,'convertfrom','datenum');
    
    % QC adcp
    aries.vel(repmat(min(aries.corr,[],2)<60,[1 size(aries.vel,2),1])) = nan;
    polly.vel(repmat(min(polly.corr,[],2)<100,[1 size(polly.vel,2) 1])) = nan;
    
    if i == 1
        polly.vel(abs(polly.vel)>2) = nan;
    end
    
    % ms depth
    [~,polly.idx_ms] = min(abs(polly.cell_depth-ms_depth(i)));
    [~,aries.idx_ms] = min(abs(aries.cell_depth-ms_depth(i)));
    
    % distance to meltstake
    polly.ms_dist = gcDist(polly.vessel_lat,polly.vessel_lon,ms_loc(1),ms_loc(2));
    aries.ms_dist = gcDist(aries.vessel_lat,aries.vessel_lon,ms_loc(1),ms_loc(2));
    
    % distance to terminus
    polly.term_dist =  nan(1,length(polly.nuc_time));
    aries.term_dist =  nan(1,length(aries.nuc_time));
    for j = 1:length(polly.term_dist)
        polly.term_dist(j) = min(gcDist(term.lat,term.lon,polly.vessel_lat(j),polly.vessel_lon(j)));
    end
    for j = 1:length(aries.term_dist)
        aries.term_dist(j) = min(gcDist(term.lat,term.lon,aries.vessel_lat(j),aries.vessel_lon(j)));
    end
    
    % meltstake proximity
    ms_prox_lim = [50*tand(25) 1200];
    polly.idx_ms_prox = polly.ms_dist>=ms_prox_lim(1) & polly.ms_dist<=ms_prox_lim(2);
    aries.idx_ms_prox = aries.ms_dist>=ms_prox_lim(1) & aries.ms_dist<=ms_prox_lim(2);
    
    % terminus proximity
    term_prox_lim = [0 200];
    polly.idx_term_prox = polly.term_dist>=term_prox_lim(1) & polly.term_dist<=term_prox_lim(2);
    aries.idx_term_prox = aries.term_dist>=term_prox_lim(1) & aries.term_dist<=term_prox_lim(2);
    
    % combine ms and terminus proximity
    polly.idx_prox = polly.idx_ms_prox & polly.idx_term_prox;
    aries.idx_prox = aries.idx_ms_prox & aries.idx_term_prox;
    
    % grid the velocities
    grid_dx = 25; % m
    grid_dlat = grid_dx/(2*pi*6.4e6/360);
    grid_dlon = grid_dlat/cosd(mean(polly.vessel_lat));
    p_loc = [polly.vessel_lat(polly.idx_prox)' polly.vessel_lon(polly.idx_prox)'];
    a_loc = [aries.vessel_lat(aries.idx_prox)' aries.vessel_lon(aries.idx_prox)'];
    p_vel = polly.vel(:,:,polly.idx_prox);
    a_vel = aries.vel(:,:,aries.idx_prox);
    lat_grid = (floor(min([p_loc(:,1);a_loc(:,1)])/grid_dlat):ceil(max([p_loc(:,1);a_loc(:,1)])/grid_dlat))'*grid_dlat;
    lon_grid = (floor(min([p_loc(:,2);a_loc(:,2)])/grid_dlon):ceil(max([p_loc(:,2);a_loc(:,2)])/grid_dlon))'*grid_dlon;
    nlat = length(lat_grid)-1;
    nlon = length(lon_grid)-1;
    polly.uv_grid = nan(nlat,nlon,2,2);
    aries.uv_grid = nan(nlat,nlon,2,2);
    polly.n_grid = nan(nlat,nlon,2);
    aries.n_grid = nan(nlat,nlon,2);
    for j = 1:nlat
        lat_lim = lat_grid(j:j+1);
        for k = 1:nlon
            lon_lim = lon_grid(k:k+1);
            p_idx = p_loc(:,1)>=lat_lim(1) & p_loc(:,1)<lat_lim(2) & p_loc(:,2)>=lon_lim(1) & p_loc(:,2)<lon_lim(2);
            a_idx = a_loc(:,1)>=lat_lim(1) & a_loc(:,1)<lat_lim(2) & a_loc(:,2)>=lon_lim(1) & a_loc(:,2)<lon_lim(2);
            polly.n_grid(j,k) = sum(p_idx);
            aries.n_grid(j,k) = sum(a_idx);
            % sfc
            polly.uv_grid(j,k,1,1) = mean(p_vel(1,1,p_idx),'all','omitnan');
            polly.uv_grid(j,k,2,1) = mean(p_vel(1,2,p_idx),'all','omitnan');
            aries.uv_grid(j,k,1,1) = mean(a_vel(1,1,a_idx),'all','omitnan');
            aries.uv_grid(j,k,2,1) = mean(a_vel(1,2,a_idx),'all','omitnan');
            % ms depth
            polly.uv_grid(j,k,1,2) = mean(p_vel(polly.idx_ms,1,p_idx),'all','omitnan');
            polly.uv_grid(j,k,2,2) = mean(p_vel(polly.idx_ms,2,p_idx),'all','omitnan');
            aries.uv_grid(j,k,1,2) = mean(a_vel(aries.idx_ms,1,a_idx),'all','omitnan');
            aries.uv_grid(j,k,2,2) = mean(a_vel(aries.idx_ms,2,a_idx),'all','omitnan');
        end
    end
    uv_grid = mean(cat(5,polly.uv_grid,aries.uv_grid),5,'omitnan');
    n_grid = polly.n_grid+aries.n_grid;
    
    % plot
    q_scale = 4*grid_dlat;
    lw = 1;
    fig = figure(1);
    ax = subplot(1,3,i);%axes(fig);
    clear h
    hold on
    box on
    
%     [LAT,LON] = meshgrid(lat_grid,lon_grid);
    % plot(LON,LAT,'+','color',0.5*[1 1 1])
    
    % terminus line
    h(1) = patch(term.lon,term.lat,0.2*[1 1 1],'facealpha',0.25,'edgecolor','k','linewidth',lw);
    % plot(term.lon,term.lat,'k-','linewidth',lw)
    
    % RHIB tracks
    % plot(polly.vessel_lon,polly.vessel_lat,'color',0.5*[1 1 1],'linewidth',lw-0.2)
    % plot(aries.vessel_lon,aries.vessel_lat,'color',0.5*[1 1 1],'linewidth',lw-0.2)
    % plot(polly.vessel_lon(polly.idx_prox),polly.vessel_lat(polly.idx_prox),'.','color',colors(1),'linewidth',lw)
    % plot(aries.vessel_lon(aries.idx_prox),aries.vessel_lat(aries.idx_prox),'.','color',colors(2),'linewidth',lw)
    
    % velocity arrows
    % quiver([polly.vessel_lon(polly.idx_prox) aries.vessel_lon(aries.idx_prox)]',[polly.vessel_lat(polly.idx_prox) aries.vessel_lat(aries.idx_prox)]',[squeeze(polly.vel(1,1,polly.idx_prox)); squeeze(aries.vel(1,1,aries.idx_prox))],[squeeze(polly.vel(1,2,polly.idx_prox)); squeeze(aries.vel(1,2,aries.idx_prox))],1.2,'r')
    h(4) = quiver(lon_grid(1:end-1)-0.*grid_dlon,lat_grid(1:end-1)-0.*grid_dlat,uv_grid(:,:,1,2)*q_scale/cosd(56.8),uv_grid(:,:,2,2)*q_scale,'off','r','linewidth',1);
    h(3) = quiver(lon_grid(1:end-1)-0.*grid_dlon,lat_grid(1:end-1)-0.*grid_dlat,uv_grid(:,:,1,1)*q_scale/cosd(56.8),uv_grid(:,:,2,1)*q_scale,'off','k','linewidth',1);
    
    % ms location
    h(2) = plot(ms_loc(2),ms_loc(1),'ro','markersize',4,'markerfacecolor','r','linewidth',lw);
    
    legend(h,{'glacier','ms location','sfc',sprintf('ms depth (%d m)',ms_depth(i))},'location','northeast','autoupdate','off')
    % quiver(-132.359,56.835,-1*q_scale/cosd(56.8),0*q_scale,'off','k','linewidth',1);
    % quiver(-132.359,56.835,-0*q_scale/cosd(56.8),1*q_scale,'off','k','linewidth',1);
    
    xlim(-132-[0.367 0.358])
    ylim([56.832 56.84])
    pbaspect([diff(xlim)*cosd(56.8)/diff(ylim) 1 1])

    % save data
    grid_vel(i).term = term;
    grid_vel(i).ms_loc = ms_loc;
    grid_vel(i).lat = lat_grid;
    grid_vel(i).lon = lon_grid;
    grid_vel(i).vel = uv_grid;
    grid_vel(i).n = n_grid;
    grid_vel(i).depth = {polly.cell_depth([1 polly.idx_ms]), aries.cell_depth([1 aries.idx_ms])};
end

% save
% save platform_comparison_data\grid_vel.mat grid_vel