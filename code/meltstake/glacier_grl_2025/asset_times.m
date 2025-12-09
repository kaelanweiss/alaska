% Script to plot out when Polly, RiffRaft, MS, etc were in the water in
% 2024
%
% KJW
% 8 Dec 2025
clear

load('G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Polly\deployments_all\adcp_combo_QC.mat')
load('G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2406\data\processed\Polly\deployments_all\old\ctd_combo.mat')

dep_tbl = loadMSInfo(26:28);

load glacier_clrs.mat

% just keep glacier data
idx_adcp = false(length(adcp),1);
for i = 1:length(adcp)
    idx_adcp(i) = adcp(i).nuc_time(1) > datenum(2024,7,1);
end
adcp = adcp(idx_adcp);

idx_ctd = false(length(ctd),1);
for i = 1:length(ctd)
    if ~isempty(ctd(i).time)
        idx_ctd(i) = ctd(i).time(1) > datenum(2024,7,1);
    end
end
ctd = ctd(idx_ctd);


%% plot
lw = 1;
fig = figure(1); clf
hold on

% MS deps
for i = 1:3
    seg_tbl = loadMSInfo(dep_tbl.Number(i),'segments');
    for j = 1:size(seg_tbl,1)
        plot([seg_tbl.Start(j) seg_tbl.End(j)],[1 1],'-','linewidth',lw,'color',dep_clrs(i,:))
    end
end

% Polly ADCP
for i = 1:length(adcp)
    plot(datetime(extrema(adcp(i).nuc_time),'convertfrom','datenum'),2*[1 1],'k-','linewidth',lw)
end

% Polly CTD
for i = 1:length(ctd)
    plot(datetime(ctd(i).time,'convertfrom','datenum'),2+0*ctd(i).time,'rx','linewidth',lw)
end


ylim([0 3])