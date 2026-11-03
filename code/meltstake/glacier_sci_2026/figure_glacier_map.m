% Script to plot terminus positions with deployment locations, approx plume
% locations, bathy, etc
%
% KJW
% 28 Jun 2026
clear

dep_tbl = loadMSInfo(26:28);
load glacier_clrs.mat

land = readgeotable("landareas.shp");

% multibeam from Eidam et al 2020
mb_file = '../../../../../code/map_plotting/data/xs_multibeam/LeConte_2018_seafloorbathymetry.nc';
mb = struct();
for fld = {'x','y','z'}
    mb.(fld{1}) = ncread(mb_file,fld{1});
end
% trim
mb.z(mb.z < -5000) = nan;
x_idx = any(~isnan(mb.z'));
y_idx = any(~isnan(mb.z));
mb.x = mb.x(x_idx,y_idx);
mb.y = mb.y(x_idx,y_idx);
mb.z = mb.z(x_idx,y_idx);

% terminus files
term_files = {'G:Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Drone\LeconteTerminus\07122024\02\TermLatLon\LeconteTermVert07122024_2.csv',...
              'G:Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Drone\LeconteTerminus\07152024\01\TermLatLon\LeconteTermVert07152024_1.csv',...
              'G:Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Drone\LeconteTerminus\07162024\02\TermLatLon\LeconteTermVert07162024_2.csv'};
nt = length(term_files);
term(nt) = struct('lat',[],'lon',[]);
for i = 1:nt
    termi = readtable(term_files{i});
    term(i).lat = termi.Y;
    term(i).lon = termi.X;

    % convert to UTM
    [term(i).x,term(i).y] = ll2utm(termi.Y,termi.X,'wgs84',8);

    if any(i==[2 3])
        term(i).lat = hannFilter(term(i).lat,3);
        term(i).lon = hannFilter(term(i).lon,3);
    end
end

% trim mb data that's behind 2024 terminus locations
ext_pts = [662875 6304607;661760 6301438];
in_glac = false(size(mb.z));
for i = 1:nt
    in_glac = in_glac | inpolygon(mb.x,mb.y,[term(i).x; ext_pts(:,1)],[term(i).y; ext_pts(:,2)]);
end
mb.z(in_glac) = nan;

% dep locations
% ms_loc = nan(nt,2);
% for i = 1:nt
%     C = strsplit(dep_tbl.Location{i},', ');
%     ms_loc(i,1) = str2double(C{1});
%     ms_loc(i,2) = str2double(C{2});
% end
ms_loc = [dep_tbl.Latitude dep_tbl.Longitude];
ms_loc(1,:) = ms_loc(1,:)+40*[cosd(ms_loc(1,1)) -1]*360/(6400000*2*pi); % jitter
ms_loc_utm = nan(size(ms_loc));
[ms_loc_utm(:,1),ms_loc_utm(:,2)] = ll2utm(ms_loc(:,1),ms_loc(:,2),'wgs84',8);

% CTD locations
load F:/meltstake/data/proc/glacier_2024_ctd.mat

% terminus raster
% term_rast = readgeoraster('G:Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Drone\LeconteTerminus\07162024\02\DEM07162024HighRes_02.tif');

% sentinel image
[Z,R] = readgeoraster('../../../../../code/map_plotting/images/2024-07-08-00_00_2024-07-08-23_59_Sentinel-2_L2A_True_color.tiff');
[lat_sat,lon_sat] = geographicGrid(R);
[x_sat,y_sat] = ll2utm(lat_sat,lon_sat,'wgs84',8);

% alaska coast
coast = readgeotable('../../../../../code/map_plotting/data/Alaska_1_250,000/Alaska_1%3A250%2C000.shp');

%% plot
% fs = 11;
% lw = 1.5;
% marker_scale = 40;
% 
% fig = figure(1); clf
% setFigureSize(fig,[7 6]);
% 
% p = projcrs(26708);
% mx = mapaxes(fig,'projectedcrs',p,'ticklabelformat','dd','position',axgridpos(1,1,1,[0 0 .15 .15],[.1 -.03]));
% hold(mx,'on')
% set(mx,'fontsize',fs)
% set(mx,'graticulealpha',.05)
% set(mx.Scalebar)
% 
% % geoplot(mx,coast(coast.OBJECTID==7155,:))
% % geoshow(mx,lat_sat,lon_sat,Z)
% 
% % for i = 1:nt
% %     geoplot(term(i).lat,term(i).lon,'-','color',dep_clrs(i,:),'linewidth',lw)
% % end
% % for i = 1:nt
% %     ms = marker_scale;
% %     if strcmp(smbs{i},'square')
% %         ms = ms*1.3;
% %     end
% %     geoscatter(ms_loc(i,1),ms_loc(i,2),ms,dep_clrs(i,:),'filled','linewidth',lw,'marker',smbs{i},'markeredgecolor','k')
% % end
% 
% mx2 = mapaxes(fig,'projectedcrs',projcrs(102016,"Authority","ESRI"),'ticklabelformat','dd','position',axgridpos(3,2,1,[0 0 .02 .02],[.0 -.0]));
% geoplot(mx2,land,'facecolor',colors(1))

% %% satellite raster
% figure(2); clf
% % p2 = projcrs(26708)
% % [x_sat,y_sat] = projfwd(p,lat_sat,lon_sat);
% % mapshow(x_sat,y_sat,Z)
% geoshow(Z,R);
% hold on
% for i = 1:nt
%     plot(term(i).lon,term(i).lat,'-','color',dep_clrs(i,:),'linewidth',lw)
% end
% for i = 1:nt
%     ms = marker_scale;
%     if strcmp(smbs{i},'square')
%         ms = ms*1.3;
%     end
%     scatter(ms_loc(i,2),ms_loc(i,1),ms,dep_clrs(i,:),'filled','linewidth',lw,'marker',smbs{i},'markeredgecolor','k')
% end

%% janky pcolor
fs = 10;
lw = 1.2;
marker_scale = 25;

fig_scale = 1.2;
fig3 = figure(3); clf
setFigureSize(fig3,fig_scale*[5 5])
ax_pos = axgridpos(1,1,1,[0 0 .18 .18],fig_scale*[-.02 -.02]);
ax_sat = axes(fig3,'position',ax_pos);
ax_mb = axes(fig3,'position',ax_pos);

p0 = [660200+20 6301750-120];

% satellite
Z_gray = rgb2gray(double(Z)/(2^16));
pcolor(ax_sat,x_sat-p0(1),y_sat-p0(2),Z_gray)
shading(ax_sat,'flat')
colormap gray

% mb
pcolor(ax_mb,mb.x-p0(1),mb.y-p0(2),double(-mb.z))
shading(ax_mb,'flat')
% colormap(ax_mb,'parula')
cmocean('-tempo',ax_mb)
hold(ax_mb,'on')

% set mb axes in sat axes
ax_mb.Color = 'none';
linkaxes([ax_sat ax_mb])

% terminus position
for i = 1:nt
    plot(ax_mb,term(i).x-p0(1),term(i).y-p0(2),'-','color',dep_clrs(i,:),'linewidth',1.2*lw)
end

% deployment locations
for i = 1:nt
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    scatter(ax_mb,ms_loc_utm(i,1)-p0(1),ms_loc_utm(i,2)-p0(2),ms,dep_clrs(i,:),'filled','linewidth',0.8*lw,'marker',smbs{i},'markeredgecolor','k')
end

% limits
L = 1300;
xlim(ax_mb,[0 L])
ylim(ax_mb,[0 L])
clim(ax_mb,[-230 0])
clim(ax_sat,[0.3 1])

% ticks and labels
set(ax_sat,'xtick',[],'ytick',[])
set(ax_mb,'xtick',0:500:L,'ytick',0:500:L)
set(ax_mb,'yticklabelrotation',90)
xlabel(ax_mb,'easting [m]','fontsize',fs)
ylabel(ax_mb,'northing [m]','fontsize',fs)

% cbar
ax_mb_pos = ax_mb.Position;
cbar = colorbar(ax_mb,'orientation','horizontal');
cbar.AxisLocationMode = 'manual';
cbar.AxisLocation = 'in';
cbar.Position(1) = .19;
cbar.Position(2) = sum(ax_mb_pos([2])) + 0.03;
cbar.Position(3) = 0.22;
cbar.Position(4) = .04;
cbar.Label.String = 'depth [m]';
cbar.Label.FontSize = 8;
cbar.Color = 0.85*[1 1 1];

ctd_version = 1;
if ctd_version
    % ctd locations
    for i = 1:nt
        [x_ctd,y_ctd] = ll2utm(ctd(i).lat,ctd(i).lon,'wgs84',8);
        scatter(ax_mb,x_ctd-p0(1),y_ctd-p0(2),6,dep_clrs(i,:),'filled','markeredgecolor','k')
    end
    cbar.Position(2) = sum(ax_mb.Position([2 4])) + .02;
    cbar.AxisLocation = 'out';
    cbar.Color = 0.15*[1 1 1];
    clim(ax_mb,[-280 0])
end

% CTD text
% scatter(ax_mb,440,1040,6,'w','filled','markeredgecolor','k')
% text(ax_mb,440+30,1040,6,'CTD','fontsize',8,'fontweight','bold', ...
%     'horizontalalignment','left','verticalalignment','middle')

%%
% fs = 11;
% lw = 1.5;
% marker_scale = 40;
% 
% fig = figure(4); clf
% setFigureSize(fig,8*[1 1]);
% 
% lat0 = 58.30;
% lon0 = -134.42;
% 
% p = projcrs(102037,'authority','esri');
% p.ProjectionParameters.LatitudeOfNaturalOrigin = lat0;
% p.ProjectionParameters.LongitudeOfNaturalOrigin = lon0;
% 
% mx2 = mapaxes(fig,'projectedcrs',p);
% geoplot(mx2,land,'facecolor',0.5*[1 1 1])

%%
% figure(5)
% axesm("ortho", ...
%     "Origin", [58.30 -134.42 0], ...
%     "Frame", "on", ...
%     "Grid", "on");
% 
% axis off
% land = shaperead("landareas.shp", "UseGeoCoords", true);
% 
% geoshow(land, ...
%     "FaceColor", [0.85 0.88 0.80], ...
%     "EdgeColor", [0.35 0.35 0.35]);
% 
% plotm(58.30, -134.42, "r.", "MarkerSize", 20);