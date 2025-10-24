% Script to manually define the location of the ice edge in ADCP beams for
% each segment. These higher-resolution ice edge estimates are used for
% remapping beam velocities into an ice-aware reference frame.
%
% KJW
% 23 Oct 2025

clear

raw_dir = 'F:meltstake/data/raw';
proc_dir = 'F:meltstake/data/proc';

% choose deployment
dep_num = 26;
seg_tbl = loadMSInfo(dep_num,'segments');
ns = size(seg_tbl,1);

% load data
load(fullfile(raw_dir,seg_tbl.Folder{1},'adcp','adcp.mat'))
adcp.burst.range = cast(adcp.burst.range,'double'); % make sure adcp.burst.range is double precision >:(
adcp.burst.time = datenum(adcp.burst.time);
idxr = adcp.burst.range>=0.3 & adcp.burst.range<=(max(seg_tbl.rmax)+0.3);

%% set up figure
fig = figure(1); clf
clear ax
pad = [0 0.02 0.085 0.07];
shift = [-0.02 0];
for i = 1:adcp.burst.nbeams
    ax(i) = axes(fig,'position',axgridpos(adcp.burst.nbeams,1,i,pad,shift));
    pcolor(adcp.burst.time,adcp.burst.range(idxr),adcp.burst.vel(:,idxr,i)')
    shading flat
    cmocean('bal',ax(i))
    clim(ax(i),.08*[-1 1])
end
linkaxes(ax)

% each beam in each segment needs a matrix with (relative) timestamps and
% ranges
adcp_ice_edges = cell(ns,1);
for i = 1:ns
    adcp_ice_edges{i} = cell(adcp.burst.nbeams,1);
end

% loop through segments, clicking location of furthest ocean adcp cell
for i = 1:ns
    fprintf('Segment %d/%d\n',i,ns)
    xlim(ax(1),datenum([seg_tbl.Start(i) seg_tbl.End(i)]))
    datetick(ax(3),'keeplimits')

    for j = 1:adcp.burst.nbeams
        fprintf('\tbeam %d\n',j)
        [tg,rg] = ginput();
        cell_nums = ceil((rg-adcp.burst.range(1))/adcp.burst.cellsize)+1;
        r_ice = adcp.burst.range(cell_nums)';
        adcp_ice_edges{i}{j} = [tg r_ice];
    end
end

% save?
save_input = input('save file (y/n)? ','s');
if strcmpi(save_input,'y')
    save(fullfile(proc_dir,seg_tbl.Folder{1},'adcp_ice_edges.mat'),'adcp_ice_edges')
end
