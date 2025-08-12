% Script to do a quick inspection of the meltstake deployments from 2024

clear

raw_path = 'F:AK_202406/data/meltstake';

d = dir(fullfile(raw_path,'ms*'));

% remove folders that aren't successful deployments
good_deps = false(length(d),1);
dt_start = NaT(size(good_deps));
for i = 1:length(d)
    good_deps(i) = ~isnan(str2double(d(i).name(end-3:end)));
    if good_deps(i)
        dt_start(i) = datetime(d(i).name(6:end),'format','yyyyMMdd_HHmm');
    end
end
d = d(good_deps);
dt_start = dt_start(good_deps);
ndeps = sum(good_deps);

% sort by date
[dt_start,sidx] = sort(dt_start);
d = d(sidx);

% load and plot (ADCP)
for i = 1:ndeps
    fprintf('=== %s ===\n',d(i).name)
    load(fullfile(raw_path,d(i).name,'adcp','adcp.mat'))

    % plot
    dep_ttl = strrep(d(i).name,'_','\_');
    ax1 = adcpQuickPlot(figure(2*i-1),adcp,'vel',0.15,NaT,[0 1],8);
    title(ax1(1),dep_ttl)

    figure(2*i); clf
    ax2 = axes(figure(2*i));
    plot(adcp.attitude.time,adcp.attitude.pressure,'k-');
    title([dep_ttl ' pressure'])

    linkaxes([ax1 ax2],'x')
end