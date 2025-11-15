% Script to take a quick look at the calibrated T-rake data from 2024
%
% KJW
% 29 Oct 2025
clear

raw_dir = 'F:meltstake/data/raw';
load glacier_grl_2025\glacier_clrs.mat

% deployment info
dep_tbl = loadMSInfo(28);
dep_nums = dep_tbl.Number;
ndeps = length(dep_nums);

% plotting colors
load glacier_grl_2025\glacier_clrs.mat
cmap = flipud(cmocean('hal',16));

% loop through deployments
for i = 1:ndeps
    % load data
    trake = loadTrake(dep_nums(i),dep_tbl);
    load(fullfile(raw_dir,dep_tbl.Folder{i},'rbr','T.mat'))

    % create figure
    fig = figure(1); clf
    ax = axes(fig);
    hold(ax,'on'); box(ax,'on')
    clear h
    h_lbls = cell(length(T)+length(trake.good_channels),1);
    % plot solos
    for j = 1:length(T)
        idxt = T(j).time>=dep_tbl.Start(i) & T(j).time<=dep_tbl.End(i);
        h(length(T)-j+1) = plot(ax,T(j).time(idxt),T(j).values(idxt),'-','color',T_clrs(length(T)-j+1,:),'linewidth',1);
        h_lbls{length(T)-j+1} = sprintf('%d mm',100+150*(j-1));
    end

    % plot Trake
    for j = 1:length(trake.good_channels)
        chan = trake.good_channels(j);
        colorj = cmap(floor(size(cmap,1)/size(trake.T,2))*chan,:);
        h(length(T)+length(trake.good_channels)-j+1) = plot(ax,trake.time(1:10:end),trake.T(1:10:end,chan),'-','color',colorj);
        h_lbls{length(T)+length(trake.good_channels)-j+1} = sprintf('%d mm',trake.dist(chan));
    end
    
    % legend
    lgd = legend(ax,h,h_lbls,'location','eastoutside','fontsize',9);

    % labels
    title(strrep(dep_tbl.Folder{i},'_','\_'))
    ylabel(ax,'temperature [^\circC]')
end




% % % % %
% Function to read out time, dist, and T from Trake file
function trake = loadTrake(dep_num,dep_tbl)
    trake = struct();
    switch dep_num
        case 26
            trake_dep = 2;
        case 27
            trake_dep = 5;
        case 28
            trake_dep = 6;
    end

    t1 = dep_tbl.Start(dep_tbl.Number==dep_num);
    t2 = dep_tbl.End(dep_tbl.Number==dep_num);

    % load
    proc_dir = 'G:Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Meltstake\Trake';
    load(fullfile(proc_dir,sprintf('trake_calibrated_deploy_%d.mat',trake_dep)),'trake_calibrated','trake_data','c','d')
    
    % save relevant data
    idxt = trake_data.data_time>=t1 & trake_data.data_time<=t2;
    trake.time = trake_data.data_time(idxt);
    trake.dist = c.distances;
    trake.T = trake_calibrated(idxt,:);
    trake.good_channels = c.good_channels+1;   
    
end
