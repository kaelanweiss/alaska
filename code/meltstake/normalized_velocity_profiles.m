% Script to plot velocity from ADV within a deployment normalized by the
% ADCP "outer" velocity.
%
% notes: probably want to process adv into segments at some point and then
% update
%
% KJW
% 22 May 2026
clear

proc_dir = 'F:/meltstake/data/proc';
raw_dir = 'F:/meltstake/data/raw';

% deployment number
dep_num = 27;
ms_tbl = loadMSInfo(dep_num,'segments');
nsegs = size(ms_tbl,1);
dep_name = ms_tbl.Folder{1};

% load data
load(fullfile(proc_dir,dep_name,'adv','pck_edges.mat'))
load(fullfile(proc_dir,dep_name,'adv','svol_in_ocean.mat'))
load(fullfile(raw_dir,dep_name,'adv','adv.mat'),'adv')
load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'))

adcp = msADCPTransform(adcp,adcp.burst.processing.cor_min,adcp.burst.processing.amp_min);
adv = msADVTransform(adv,adcp.attitude);
adv.vel_ice = [adv.vel_ice vecnorm(adv.vel_ice,2,2)];
adv = burstAverageADV(adv);
b_dt = adv.samples_per_burst/adv.sample_rate;

idxr = adcp.burst.range <= 0.2;

%% calculate velocity things
% u_outer = nan(nsegs,4,2); % segment, (u,v,w,U), (mean, stdev)

i = 10;
%for i = 1:nsegs
% load adcp
% load(fullfile(proc_dir,dep_name,sprintf('adcp%d.mat',i)))
t1 = ms_tbl.Start(i);
t2 = ms_tbl.End(i);
idxb = adv.time_avg >= t1 & adv.time_avg <= t2;
idxy = edges.time >= t1 & edges.time <= t2;
u_adv = adv.vel_avg(idxb,:);
y_adv = (edges.boundary_location(idxb)-157)/10;
nb = sum(idxb);
u_outer = nan(nb,4,2);
ti = adv.time_avg(idxb);
for j = 1:nb
    idxj = adcp.burst.time >= ti(j) & adcp.burst.time <= ti(j)+seconds(b_dt);
    velj = squeeze(mean(adcp.burst.vel_ice(idxj,idxr,[1 3 2]),2,'omitnan'));
    velj = [velj vecnorm(velj,2,2)];
    u_outer(j,:,1) = mean(velj,'omitnan');
    u_outer(j,:,2) = std(velj,'omitnan');
end

too_slow = [abs(u_outer(:,1:3))<0.01 abs(u_outer(:,4))<0.03];

u_norm = u_adv./cat(3,u_outer(:,:,1),abs(u_outer(:,:,2)));
u_norm(repmat(too_slow,[1 1 2])) = nan;

% figure(1); clf
% clear ax
comp_lbls = {'u','v','w','|u|'};
for j = 4:-1:1
    % ax(j) = subplot(4,1,j); hold on
    errorbar(ax(j),y_adv,u_norm(:,j,1),u_norm(:,j,2)/10,'ko','linewidth',0.8,'markersize',3,'markerfacecolor','k')
    % ylabel(ax(j),sprintf('%s_{adv}/%s_{adcp}',comp_lbls{j},comp_lbls{j}))
    % xlabel(ax(j),'y [cm]')
end
% linkaxes(ax,'x')
% xlim(ax(1),[-0.5 max(ceil(2*y_adv))/2])