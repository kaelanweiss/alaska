% Script to perform bin mapping of the ADCP data relative to the ice
% locations (determined using script "manual_adcp_ice_edges.m").
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
load(fullfile(proc_dir,seg_tbl.Folder{1},'adcp_ice_edges.mat'))

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