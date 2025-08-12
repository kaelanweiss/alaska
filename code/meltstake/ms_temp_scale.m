% script to calculate temp scales for meltstake deployments
% 
% NOT USED FOR ACTUAL T SCALE CALCULATION, ACTUAL SCRIPT IS
% "ms_scales_osm.m" - 30 Sep 2024

clear

% load deployment info
i = 2;

tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
mstbl = readtable(tbl_path);
mstbl = mstbl(~isnan(mstbl.Number),:);

% load rbr temperature
load(fullfile('F:meltstake\data\raw',mstbl.Folder{i},'rbr','T.mat'))

%%
ns = length(T);
T_all = nan(length(T(1).time),ns);
hann_width = 21;
for i = 1:ns
    T_all(:,i) = hannFilter(T(i).values,hann_width);
end
time = T(1).time;

% convert to datetime
if ~isdatetime(time)
    time = datetime(time,'convertfrom','datenum');
end

% figure out window length for half-overlapping windows
T_total = time(end)-time(1);
T1 = minutes(20)/5;
T2 = minutes(30)/5;

[L,n] = getWindowLength(T_total,T1,T2);

% calculate T scale
[T_scale,T_std,t_scale] = msTempScale(T_all,time,n,L);

%% plot
figure(6); clf; hold on

% original record
for i = 1:ns
    plot(time,T_all(:,i),'color',[1 1 1]*0.4 + 0.1*i)
end

% scale with std
dt_plt = L/20;
for i = 1:ns
    errorbar(t_scale+dt_plt*(i-1),T_scale(:,i),T_std(:,i),'.','color',colors(i),'linewidth',1)
    plot(t_scale+dt_plt*(i-1),T_scale(:,i),'.','markersize',12,'color',colors(i))
end

i = ns+1;
errorbar(t_scale+dt_plt*(i-1),T_scale(:,i),T_std(:,i),'k.','linewidth',1)
plot(t_scale+dt_plt*(i-1),T_scale(:,i),'k.','markersize',12)




