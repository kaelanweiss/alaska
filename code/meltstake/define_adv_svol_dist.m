% Script to manually determine the location of the ice boundary in the ADV
% backscatter data. This is used for determining the y-coordinate of
% measurements, and relies on having the ice edges detected already (in
% pck_edges.mat files in processed directory).
%
% KJW
% 21 May 2026
clear

raw_dir = 'F:/meltstake/data/raw';
proc_dir = 'F:/meltstake/data/proc';
load F:/adv/mean_profile.mat

% save flag
save_output = 1;

% deployment
dep_num = 27;

% segment info
ms_tbl = loadMSInfo(dep_num,'segments');
nsegs = size(ms_tbl,1);
dep_name = ms_tbl.Folder{1};

% load data
load(fullfile(raw_dir,dep_name,'adv','adv.mat'))
load(fullfile(proc_dir,dep_name,'adv','pck_edges.mat'))
if ~isfield(edges,'boundary_location')
    edges.boundary_location = nan(size(edges.time));
end

%% make pck amplitude
pck_amp = mean(adv.pck_amp,3);
diff_t = seconds(diff(adv.pck_time([1:3 end-2:end])));
if diff_t(1) > diff_t(2) % pad the beginning
    pck_amp = cat(1,pck_amp(1,:),pck_amp);
end
if diff_t(end-1) < diff_t(end) % pad the end
    pck_amp = cat(1,pck_amp,pck_amp(end,:));
end
pck_amp = (pck_amp(1:2:end,:) + pck_amp(2:2:end,:))/2;
pck_amp = pck_amp./(repmat(0.9*mean_prof',[size(pck_amp,1) 1])+0);

%% initial plot
% set up figure
fig = figure(1); clf
ax = axes(fig); hold on

% pcolor
pcolor(adv.burst_time,adv.pck_dist(1,:),pck_amp')
shading flat
clim([0.9 1.6])

% edges
mrks = {'o','^','square'};
for i = 1:3
    scatter(edges.time,edges.pos(:,i),8,'marker',mrks{i},'markeredgecolor','k','markerfacecolor','k')
end

% boundary if it exists
p_boundary = plot(edges.time,edges.boundary_location,'r-');

%% loop through segments
for i = 1:nsegs
    fprintf('segment %d\n',i)
    % time range
    t1 = ms_tbl.Start(i);
    t2 = ms_tbl.End(i);
    idxi = edges.time >= t1 & edges.time < t2;
    xq = edges.time(idxi);
    % adjust limits
    xlim(ax,[t1 t2]+seconds(45*[-1 1]))
    ylim(ax,mean(edges.pos(idxi,:),'all')+70*[-1 1])

    % boundary loop
    while true
        % turn off line
        p_boundary.YData(:) = nan;
        % new boundary
        fprintf('point and click boundary...\n')
        [xg,yg] = ginput();
        [xg,sidx] = sort(xg);
        yg = yg(sidx);
        y_boundary = interp1(num2ruler(xg,ax.XAxis),yg,edges.time(idxi));
        % plot
        p_boundary.XData = edges.time(idxi);
        p_boundary.YData = y_boundary;
        % check
        in = input('Continue with current boundary? [y/n] ','s');
        if strcmpi(in,'y')
            % save
            edges.boundary_location(idxi) = y_boundary;
            break
        end
    end
end

p_boundary.XData = edges.time;
p_boundary.YData = edges.boundary_location;
xlim(ax,[ms_tbl.Start(1) ms_tbl.End(end)])
ylim(ax,adv.pck_dist(1,[1 end]))

%% save output to edges file
if save_output
    save(fullfile(proc_dir,dep_name,'adv','pck_edges.mat'),'edges')
end




