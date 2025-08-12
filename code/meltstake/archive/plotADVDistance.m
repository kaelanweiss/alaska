% Script to generate figures of ADV distance measurements
clear

data_path = 'F:/AK_202309/data/meltstake';
d = dir(fullfile(data_path,'ms*'));

save_figs = false;
fig_path = 'F:/AK_202309/figures/meltstake';

ndeps = length(d);

% manual time limits
tlimits = {'2000-01-01 00:00:00','2100-01-01 00:00:00';...
           '2000-01-01 00:00:00','2100-01-01 00:00:00';...
           '2000-01-01 00:00:00','2100-01-01 00:00:00';...
           '2023-09-20 20:20:00','2023-09-20 21:10:00';...
           '2023-09-21 19:25:00','2023-09-21 20:05:00';...
           '2023-09-22 00:20:00','2023-09-22 02:00:00';...
           '','';...
           '2023-09-23 19:55:00','2023-09-24 00:05:00';...
           '2023-09-24 23:49:00','2023-09-25 00:25:00';...
           '',''};
tlimits = datetime(tlimits);

% load deployment, plot time series of distances
fpos = [-1400 -60 700*.935 700];
for i = [4 5 6 8 9]%1:ndeps
    % load if file exists
    adv_file = fullfile(d(i).folder,d(i).name,'adv','adv.mat');
    if ~exist(adv_file,'file')
        fprintf('adv file %s does not exist\n',adv_file)
        continue
    end
    fprintf('loading file %s\n',adv_file)
    load(adv_file,'adv');

    % time indexing
    hdr_idx = adv.hdr_time>=tlimits(i,1) & adv.hdr_time<=tlimits(i,2);
    hdr_time = adv.hdr_time(hdr_idx);

    dt_pressure = 16;
    pressure_idx = false(size(adv.time));
    pressure_idx(1:dt_pressure:end) = true;
    pressure_idx = pressure_idx & adv.time>=tlimits(i,1) & adv.time<=tlimits(i,2);
    depth = mean(adv.pressure(pressure_idx));

    % deployment label
    dep_lbl = sprintf('%s (depth=%.1fm)',strrep(d(i).name,'_','\_'),depth);

    % plot calculated distance
    dist = adv.dsvol_ave(hdr_idx,:)/10;
    dist(dist==0) = nan;
    figure(100+i); clf
    set(figure(100+i),'position',fpos);
    ax = gca;
    plot(hdr_time,dist,'k.-',markersize=8)
    box on
    grid on
    ylabel('dist (cm)')
    text(ax,0,1.025,sprintf('%s adv sample vol dist',dep_lbl),units='normalized')

    % plot beam distances
    figure(200+i); clf; hold on
    set(figure(200+i),'position',fpos);
    ax = gca;
    for k = 1:3
        plot(hdr_time,medianFilter(adv.dp(hdr_idx,k,1)/10,3),'.-',markersize=8,color=colors(k+3));
        plot(hdr_time,medianFilter(adv.dp(hdr_idx,k,2)/10,3),'.--',markersize=8,color=colors(k+3));
    end
    stdev = std(adv.dp(hdr_idx,:,:),1,'all')/10;
    ylim(2*[-stdev stdev]+mean(adv.dp(hdr_idx,:,:),'all')/10)
    box on
    grid on
    ylabel('dist (cm)')
    text(ax,0,1.025,sprintf('%s adv rcvr dist',dep_lbl),units='normalized')


    % plot adv pressure
%     figure(300+i); clf
%     plot(adv.time(pressure_idx),adv.pressure(pressure_idx),'k-')
%     ylabel('pressure (dbar)')

    % save figure
    if save_figs
        fig_path_dep = fullfile(fig_path,d(i).name);
        if ~exist(fig_path_dep,'dir')
            mkdir(fig_path_dep);
        end
        print(figure(100+i),fullfile(fig_path_dep,[d(i).name '_adv_dist.png']),'-dpng','-r450')
        print(figure(200+i),fullfile(fig_path_dep,[d(i).name '_adv_rcvr_dist.png']),'-dpng','-r450')
    end
end