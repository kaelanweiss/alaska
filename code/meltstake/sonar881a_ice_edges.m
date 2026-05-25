% Script to find ice edges in the raw 881a data and save output to the
% processed directory.
%
% KJW
% 21 Dec 2025
clear

raw_dir = 'F:/meltstake/data/raw';
proc_dir = 'F:/meltstake/data/proc';

% load
% dep_tbl = loadMSInfo([21:24 26:28]);
dep_tbl = loadMSInfo(26);

ndeps = size(dep_tbl,1);
sonar(ndeps) = struct('time',[],'angle',[],'range',[],'scan1',[],'scan2',[]);
aperture = 120;
for i = ndeps:-1:1
    load_sonar = load(fullfile(raw_dir,dep_tbl.Folder{i},'sonar881a','sonar881a.mat'));
    
    % trim to deployment time and forward angles
    t1 = dep_tbl.Start(i) - minutes(10);
    t2 = dep_tbl.End(i) + minutes(25);
    idxt = load_sonar.sonar.time>=t1 & load_sonar.sonar.time<=t2;
    idxp = abs(load_sonar.sonar.angle) <= aperture/2;
    load_sonar.sonar.time = load_sonar.sonar.time(idxt);
    load_sonar.sonar.angle = load_sonar.sonar.angle(idxp);
    load_sonar.sonar.scan1 = load_sonar.sonar.scan1(idxt,idxp,:);
    load_sonar.sonar.scan2 = load_sonar.sonar.scan2(idxt,idxp,:);
    
    sonar(i) = load_sonar.sonar;
end
% add some useful fields to structures
for i = 1:ndeps
    sonar(i).scan = (sonar(i).scan1+sonar(i).scan2)/2;
    sonar(i).folder = dep_tbl.Folder{i};
    [R,PHI] = meshgrid(sonar(i).range,sonar(i).angle*pi/180);
    sonar(i).X = R.*cos(PHI);
    sonar(i).Y = R.*sin(PHI);
    sonar(i).processing = struct('aperture',aperture);
end
clear load_sonar

%% do some averaging/smoothing
r_hann = [[30 30 30 30] [30 30 50]]; % radial Hanning filter
a_hann = [2*[5 5 5 5] 2*[5 7 9]]; % angular Hanning filter
tic
for i = 1:ndeps
    sonar(i).scan_smth = sonar(i).scan;
    % get size
    [nt,np,nr] = size(sonar(i).scan);
    % loop through time
    for j = 1:nt
        % smooth across range
        for k = 1:np
            % smooth
            sonar(i).scan_smth(j,k,:) = hannFilter(squeeze(sonar(i).scan(j,k,:)),r_hann(i));
        end
        % smooth across angle
        for k = 1:nr
            sonar(i).scan_smth(j,:,k) = hannFilter(squeeze(sonar(i).scan_smth(j,:,k)),a_hann(i));
        end
    end
    % save processing info
    sonar(i).processing.r_hann = r_hann(i);
    sonar(i).processing.a_hann = a_hann(i);
end
toc

%% backwards diff operation (to get amplitude and diff peaks to align)
for i = 1:ndeps
    % take forward diff
    sonar(i).diffY = diff(sonar(i).scan_smth,[],3);
    % pad an array of zeroes at r=0 to make it backwards diff
    sonar(i).diffY = cat(3,zeros([size(sonar(i).scan,[1 2]) 1]),sonar(i).diffY);
end

%% egde-finding
% set x windows for each deployment
x_min = [[0.1 0.1 0.1 0.1] [0.3 0.44 0.3]];
x_max = [[1.2 1.2 1.2 1.2] [1.1 0.8 1.2]];
for i = 1:ndeps
    % get size
    [nt,np,~] = size(sonar(i).scan);
    % create x mask
    x_mask = sonar(i).X>=x_min(i) & sonar(i).X<=x_max(i);
    sonar(i).processing.x_range = [x_min(i) x_max(i)];
    % define edge-finding field
    sonar(i).E = sonar(i).diffY.*sonar(i).scan_smth;
    sonar(i).E(~permute(repmat(x_mask,[1 1 nt]),[3 1 2])) = 0;
    sonar(i).processing.desc = 'maximize (forward difference * amplitude)';
    % find max(E) position along range dimension
    [~,imax] = max(sonar(i).E,[],3);
    sonar(i).idx_edge = imax;
    % find range values
    sonar(i).r_edge = sonar(i).range(sonar(i).idx_edge);
end

% extract boundary positions
for i = 1:ndeps
    [nt,np,~] = size(sonar(i).scan);
    % calculate positions
    sonar(i).x_edge = sonar(i).r_edge.*cosd(repmat(sonar(i).angle',[nt 1]));
    sonar(i).y_edge = sonar(i).r_edge.*sind(repmat(sonar(i).angle',[nt 1]));
end

%% quick plot
fig1 = figure(1); clf
clear ax1
% set up axes
pad1 = [.05 .05 .1 .1];
for i = 1:ndeps
    ax1(2*i-1) = axes(fig1,'position',axgridpos(2,ndeps,2*i-1,pad1,[0 0],'flip'));
    ax1(2*i) = axes(fig1,'position',axgridpos(2,ndeps,2*i,pad1,[0 0],'flip'));
    hold(ax1(2*i-1),'on')
    hold(ax1(2*i),'on')
end

% plot
i_scan = 15;
% for i_scan = 42:100
for i = 1:ndeps
    ii = min(i_scan,length(sonar(i).time));

    cla(ax1(2*i-1)); cla(ax1(2*i));
    % pcolors
    s1 = squeeze(sonar(i).scan(ii,:,:));
    % s2 = squeeze(sonar(i).scan_smth(i_scan,:,:));
    % s2 = squeeze(sonar(i).diffY(i_scan,:,:));
    s2 = squeeze(sonar(i).diffY(ii,:,:)).*s1;
    % s2 = squeeze(sonar(i).E(ii,:,:));
    pcolor(ax1(2*i-1),sonar(i).X,sonar(i).Y,s1);
    pcolor(ax1(2*i),sonar(i).X,sonar(i).Y,s2);

    x_edge = sonar(i).r_edge(ii,:)'.*cosd(sonar(i).angle);
    y_edge = sonar(i).r_edge(ii,:)'.*sind(sonar(i).angle);
    plot(ax1(2*i-1),x_edge,y_edge,'k.')
    plot(ax1(2*i),x_edge,y_edge,'k.')
    title(ax1(2*i-1),sprintf('scan %d',ii))

    for j = [2*i-1 2*i]
        shading(ax1(j),'flat')
        axis(ax1(j),'equal')
        if ~mod(j,2)
            cmocean('bal',ax1(j))
        end
    end
end
% input('')
% end


    
%% create and save output structures
save_output = 0;

if save_output
    for i = 1:ndeps
        sonar881a_edges = rmfield(sonar(i),{'scan1','scan2','scan','X','Y','diffY','E'});
        save(fullfile(proc_dir,dep_tbl.Folder{i},'sonar881a_edges.mat'),'sonar881a_edges')
    end
end