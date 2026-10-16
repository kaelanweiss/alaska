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
mb_file = '../../../../../code/map_plotting/data/LeConte_2018_seafloorbathymetry.nc';
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

    if any(i==[2 3])
        term(i).lat = hannFilter(term(i).lat,3);
        term(i).lon = hannFilter(term(i).lon,3);
    end
end

% dep locations
% ms_loc = nan(nt,2);
% for i = 1:nt
%     C = strsplit(dep_tbl.Location{i},', ');
%     ms_loc(i,1) = str2double(C{1});
%     ms_loc(i,2) = str2double(C{2});
% end
ms_loc = [dep_tbl.Latitude dep_tbl.Longitude];
ms_loc(1,:) = ms_loc(1,:)+12*[cosd(ms_loc(1,1)) -1]*360/(6400000*2*pi); % jitter

% plume locations


% terminus raster
% term_rast = readgeoraster('G:Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Drone\LeconteTerminus\07162024\02\DEM07162024HighRes_02.tif');

% sentinel image
[Z,R] = readgeoraster('../../../../../code/map_plotting/images/2024-07-08-00_00_2024-07-08-23_59_Sentinel-2_L2A_True_color.tiff');
[lat_sat,lon_sat] = geographicGrid(R);

% alaska coast
coast = readgeotable('../../../../../code/map_plotting/data/Alaska_1_250,000/Alaska_1%3A250%2C000.shp');


%% plot
fs = 11;
lw = 1.5;
marker_scale = 40;

fig = figure(1); clf
setFigureSize(fig,[7 6]);

p = projcrs(26708);
mx = mapaxes(fig,'projectedcrs',p,'ticklabelformat','dd','position',axgridpos(1,1,1,[0 0 .15 .15],[.1 -.03]));
hold(mx,'on')
set(mx,'fontsize',fs)
set(mx,'graticulealpha',.05)
set(mx.Scalebar)

% geoplot(mx,coast(coast.OBJECTID==7155,:))
% geoshow(mx,lat_sat,lon_sat,Z)

for i = 1:nt
    geoplot(term(i).lat,term(i).lon,'-','color',dep_clrs(i,:),'linewidth',lw)
end
for i = 1:nt
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    geoscatter(ms_loc(i,1),ms_loc(i,2),ms,dep_clrs(i,:),'filled','linewidth',lw,'marker',smbs{i},'markeredgecolor','k')
end

% mx2 = mapaxes(fig,'projectedcrs',projcrs(102016,"Authority","ESRI"),'ticklabelformat','dd','position',axgridpos(3,2,1,[0 0 .02 .02],[.0 -.0]));
% geoplot(mx2,land,'facecolor',colors(1))

%% staellite raster
figure(2); clf
% p2 = projcrs(26708)
% [x_sat,y_sat] = projfwd(p,lat_sat,lon_sat);
% mapshow(x_sat,y_sat,Z)
geoshow(Z,R);
hold on
for i = 1:nt
    plot(term(i).lon,term(i).lat,'-','color',dep_clrs(i,:),'linewidth',lw)
end
for i = 1:nt
    ms = marker_scale;
    if strcmp(smbs{i},'square')
        ms = ms*1.3;
    end
    scatter(ms_loc(i,2),ms_loc(i,1),ms,dep_clrs(i,:),'filled','linewidth',lw,'marker',smbs{i},'markeredgecolor','k')
end