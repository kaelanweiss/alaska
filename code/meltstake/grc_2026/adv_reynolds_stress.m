% Script to calculate Reynolds stresses from a full Meltstake deployment.
% This is a first crack at looking at near-boundary momentum fluxes.
%
%
% KJW
% 2 Jun 2026
clear

proc_dir = 'F:/meltstake/data/proc';
raw_dir = 'F:/meltstake/data/raw';

% deployment number
dep_num = 28;
ms_tbl = loadMSInfo(dep_num,'segments');
nsegs = size(ms_tbl,1);
dep_name = ms_tbl.Folder{1};

% load data
load(fullfile(proc_dir,dep_name,'adv','pck_edges.mat'))
load(fullfile(proc_dir,dep_name,'adv','svol_in_ocean.mat'))
load(fullfile(raw_dir,dep_name,'adv','adv.mat'),'adv')
load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'))
int_vel = load(fullfile(proc_dir,dep_name,'velocity_lowpass.mat'));

%% transform velocities (may be redundant)
switch dep_num
    case 26
        adv_roll_offset = 0;
    case 27
        adv_roll_offset = -21;
    case 28
        adv_roll_offset = 12.5;
end
adcp = msADCPTransform(adcp,adcp.burst.processing.cor_min,adcp.burst.processing.amp_min);
attitude = adcp.attitude;
attitude.roll_corrected = attitude.roll_corrected - adv_roll_offset;
adv = msADVTransform(adv,attitude);

% get rid of velocities in the ice
adv.vel_ice(~idx_ocean,:) = nan;

% nan out velocities not included in segments
idx_segs = false(size(adv.time));
for i = 1:nsegs
    t1 = ms_tbl.Start(i);
    t2 = ms_tbl.End(i);
    idx_segs = idx_segs | (adv.time>=t1 & adv.time<=t2);
end
adv.vel_ice(~idx_segs,:) = nan;

%% wall-normal coordinate adjustments based on flow direction
vel_seg_mean = nan(nsegs,6);
for i = 1:nsegs
    t1 = ms_tbl.Start(i);
    t2 = ms_tbl.End(i);
    % adcp
    idxt = adcp.burst.time>=t1 & adcp.burst.time<=t2;
    vel_seg_mean(i,1:3) = squeeze(mean(mean(adcp.burst.vel_ice(idxt,1:4,[1 3 2]),1,'omitnan'),2));
    % adv
    idxt = adv.time>=t1 & adv.time<=t2;
    vel_seg_mean(i,4:6) = mean(adv.vel_ice(idxt,:),'omitnan');
    % rotate about +z
    phii = atan2(vel_seg_mean(i,5),vel_seg_mean(i,4));
    if dep_num == 27% | dep_num == 28
        phii = 0;
    end
    adv_uv = adv.vel_ice(idxt,1) + 1i*adv.vel_ice(idxt,2);
    adv_uv2 = adv_uv*exp(-1i*phii);
    adv.vel_ice(idxt,1) = real(adv_uv2);
    adv.vel_ice(idxt,2) = imag(adv_uv2);
end
phi_adcp = atan2(vel_seg_mean(:,2),vel_seg_mean(:,1));
phi_adv = atan2(vel_seg_mean(:,5),vel_seg_mean(:,4));

% fudge 5 deg rotation for near-wall deployments, shouldn't affect further
% ones much
% adv_uv = adv.vel_ice(:,1)+1i*adv.vel(:,2);
% adv_uv = adv_uv/exp(1i*(-0)*pi/180);
% adv.vel_ice(:,1) = real(adv_uv);
% adv.vel_ice(:,2) = imag(adv_uv);

%% interpolate outer velocity onto adv time axis
tau_idx = 5;
tau = int_vel.tau(tau_idx); fprintf('Averaging timescale: %.1f s\n',tau)
adv.vel_outer = interp1(int_vel.time,int_vel.vel_lowpass(:,1,tau_idx+1)+1i*int_vel.vel_lowpass(:,3,tau_idx+1),adv.time);
adv.vel_angle = angle(adv.vel_outer);

% select vel range
vel_range = [.08 .12];
idx_vel_range = abs(adv.vel_outer)>vel_range(1) & abs(adv.vel_outer)<vel_range(2);
adv.vel_ice(~idx_vel_range,:) = nan;

% unwind velocity
adv_uw = adv.vel_ice(:,1) + 1i*adv.vel_ice(:,3);
adv_rot = adv_uw./exp(1i*adv.vel_angle);
adv.vel_rot = [real(adv_rot) adv.vel_ice(:,2) imag(adv_rot)];

%% Reynolds decomposition
% get ADV data on a continuous time vector
t = (adv.time(1):seconds(1/adv.sample_rate):adv.time(end))';
vel = nan(size(t,1),3);
burst_nums = unique(adv.burst_num);
orig_data = false(size(t));
for i = 1:length(burst_nums)
    idxi = adv.burst_num == burst_nums(i);
    ns = sum(idxi);
    t0i = adv.time(find(idxi,1));
    idx0 = find(t==t0i);
    if length(idx0) ~= 1
        warning('iter %d has %d starting indices',i,length(idx0))
    end
    vel(idx0:idx0+ns-1,:) = adv.vel_rot(idxi,:);
    orig_data(idx0:idx0+ns-1) = true;
end

% just a tiny bit of QC
vel(abs(vel)>0.8) = nan;

% smooth at a timescale and decompose
vel_bar = vel;
tau_rd = tau/1.;
for i = 1:3
    vel_bar(:,i) = hannFilter(vel(:,i),ceil(adv.sample_rate*tau_rd));
end
vel_prm = vel(orig_data,:) - vel_bar(orig_data,:);

% correlations
corr = nan(size(vel_prm,1),3,3);
for i = 1:3
    for j = 1:3
        corr(:,i,j) = vel_prm(:,i).*vel_prm(:,j);
    end
end
adv.corr = [corr(:,2,1) corr(:,3,1) corr(:,3,2) corr(:,1,1) corr(:,2,2) corr(:,3,3)];
adv.tke = corr(:,1,1)+corr(:,2,2)+corr(:,3,3);

% average by burst
adv.vel_all_flds = {'outer_adcp','u_rot','v_rot','w_rot','uv','uw','vw','uu','vv','ww','tke','outer_mag'};
adv.vel_all = [adv.vel_outer adv.vel_rot adv.corr adv.tke abs(adv.vel_outer)];
adv = burstAverageADV(adv,'vel_all');
adv.y_avg = [nan; (edges.boundary_location-157)/10]; % this is not a good line of code

% normalized values
vel_mag = adv.vel_avg(:,12);
vel_mag(vel_mag<prctile(vel_mag,5)) = nan;
u_norm = adv.vel_avg(:,2:4)./vel_mag;
tke = adv.vel_avg(:,11);
tke(tke<prctile(tke,5)) = nan;
uv_norm = adv.vel_avg(:,5:7)./tke;
% uv_norm = adv.vel_avg(:,5:7)./abs(adv.vel_avg(:,1)).^2;
Cd = -(adv.vel_avg(:,5)+0*adv.vel_avg(:,7))./abs(adv.vel_avg(:,12)).^2;
Cd(abs(adv.vel_avg(:,1))<.0 & adv.y_avg>500) = nan;
Cd(abs(Cd)>.1) = nan;
Cd_all = Cd;

% bin in distance
[y_adv,sidx] = sort(adv.y_avg);
sidx(isnan(y_adv)) = [];
y_adv(isnan(y_adv)) = [];
[u_bin,y_bin] = time_bin(y_adv,[adv.vel_avg(sidx,1) u_norm(sidx,:) uv_norm(sidx,:) adv.vel_avg(sidx,8:11)]',0.5,'omitnan');
u_bin = u_bin';

uv_sgf = sigmaFilter(adv.vel_avg(:,5),1.5,1,1);
[u_bin2,y_bin2] = time_bin(y_adv,[adv.vel_avg(sidx,2) uv_sgf(sidx)]',1,'omitnan');
u_bin2 = u_bin2';

% mask to 2-6 cm
idx26 = y_adv >= 2 & y_adv <= 6;
uv_26 = uv_sgf(idx26);
u_star_26 = sqrt(-uv_26);
Cd_26 = u_star_26.^2./adv.vel_avg(idx26,12).^2;

% trim to y>0.55
idx_close = y_bin < 0.55;
u_bin(idx_close,:) = nan;

%% plots
fs = 12;
figsize = [15 18];
% velocity components
vel_lbls = {{'along-flow','u/U_{ADCP}'},{'wall-normal','v/U_{ADCP}'},{'cross-flow','w/U_{ADCP}'}};
figure(11); clf
setFigureSize(figure(11),figsize);
clear ax
for i = 3:-1:1
    ax(i) = subplot(3,1,i); hold on; box on
    plot(y_bin,u_bin(:,i+1),'ko','markersize',6,'markerfacecolor',colors(1))
    plot(adv.y_avg,u_norm(:,i),'kx')
    xlabel('y [cm]','fontsize',fs)
    ylabel(vel_lbls{i},'fontsize',fs)
end
linkaxes(ax,'x')
% xlim(ax(1),[0 12])
ylim(ax(1),[-.5 2])
ylim(ax(2),0.5*[-1 1])
ylim(ax(3),0.5*[-1 1])
title(ax(1),sprintf('dep %d | T_{outer}: %.1fs | T_{Rnlds}: %.1fs | U: %.2f-%.2fm/s',dep_num,tau,tau_rd,round(vel_range(1),2),round(vel_range(2),2)),'fontsize',fs)
plbls = addPanelLabels(ax,fs-1);
for i = 1:length(plbls)
    plbls(i).String = sprintf('mean: %.2f',round(mean(u_bin(:,i+1),'omitnan'),2));
    plbls(i).Position(2) = plbls(i).Position(2)-.03;
end

% reynolds stress components
vel_lbls = {'-<u''v''>/<q>','-<u''w''>/<q>','-<v''w''>/<q>'};
figure(12); clf
setFigureSize(figure(12),figsize);
clear ax
for i = 3:-1:1
    ax(i) = subplot(3,1,i); hold on; box on
    plot(y_bin,-u_bin(:,i+4),'ko','markersize',6,'markerfacecolor',colors(2))
    plot(adv.y_avg,-uv_norm(:,i),'kx')
    xlabel('y [cm]','fontsize',fs)
    ylabel(vel_lbls{i},'fontsize',fs)
end
title(ax(1),sprintf('dep %d | T_{outer}: %.1fs | T_{Rnlds}: %.1fs | U: %.2f-%.2fm/s',dep_num,tau,tau_rd,round(vel_range(1),2),round(vel_range(2),2)),'fontsize',fs)
plbls = addPanelLabels(ax,fs-1);
for i = 1:length(plbls)
    plbls(i).String = sprintf('mean: %.4f',-round(mean(u_bin(:,i+4),'omitnan'),4));
    plbls(i).Position(2) = plbls(i).Position(2)-.03;
end

% TKE
vel_lbls = {'<q> [m^2/s^2]','<u''^2> [m^2/s^2]','<v''^2> [m^2/s^2]','<w''^2> [m^2/s^2]'};
figure(13); clf
setFigureSize(figure(13),figsize);
clear ax
for i = 4:-1:2
    ax(i) = subplot(4,1,i); hold on; box on
    plot(y_bin,u_bin(:,i+6),'ko','markersize',6,'markerfacecolor',colors(3))
    plot(adv.y_avg,adv.vel_avg(:,i+6),'kx')
    xlabel('y [cm]','fontsize',fs)
    ylabel(vel_lbls{i},'fontsize',fs)
end
ax(1) = subplot(4,1,1); hold on; box on
plot(y_bin,u_bin(:,11),'ko','markersize',6,'markerfacecolor',colors(3))
plot(adv.y_avg,adv.vel_avg(:,11),'kx')
xlabel('y [cm]','fontsize',fs)
ylabel(vel_lbls{1},'fontsize',fs)
title(ax(1),sprintf('dep %d | T_{outer}: %.1fs | T_{Rnlds}: %.1fs | U: %.2f-%.2fm/s',dep_num,tau,tau_rd,round(vel_range(1),2),round(vel_range(2),2)),'fontsize',fs)
plbls = addPanelLabels(ax,fs-1);
for i = 2:length(plbls)
    plbls(i).String = sprintf('mean: %.2e',round(mean(u_bin(:,i+6),'omitnan'),6));
    plbls(i).Position(2) = plbls(i).Position(2)-.03;
end
plbls(1).String = sprintf('mean: %.2e',round(mean(u_bin(:,11),'omitnan'),6));
plbls(1).Position(2) = plbls(i).Position(2)-.03;
linkaxes(ax)
ylim(ax(1),[0 prctile(adv.vel_avg(:,11),95)])

%% drag coefficient histogram
figure(21); clf
setFigureSize(figure(21),[8 6]);
ax = axes(figure(21));
histogram(ax,Cd_all*1000,51)
xlabel('C_d (x1000)','fontsize',fs)
xlim(ax,[-12 22])

%% mean downstream velocity and reynolds stress
pad = [.05 .08 .1 .1];
shift = [0.045 0.02];
figure(20); clf
setFigureSize(figure(20),[8.5 8]);
clear ax
for i = 2:-1:1
    ax(i) = axes(figure(20),'position',axgridpos(2,1,i,pad,shift),'fontsize',fs-2); hold on; box on
end
xlabel(ax(2),'wall-normal distance [cm]','fontsize',fs)

h20(1) = plot(ax(1),adv.y_avg,adv.vel_avg(:,2),'kx');
h20(2) = plot(ax(1),y_bin2,u_bin2(:,1),'ko','markersize',6,'markerfacecolor',colors(1));
ylim(ax(1),[0 max(get(ax(1),'YLim'))])
ylabel(ax(1),'<u> [m/s]','fontsize',fs)

plot(ax(2),adv.y_avg,-uv_sgf,'kx')
plot(ax(2),y_bin2,-u_bin2(:,2),'ko','markersize',6,'markerfacecolor',colors(2))
ylabel(ax(2),'-<u''v''> [m^2/s^2]','fontsize',fs)
% title(ax(1),sprintf('dep %d | T_{outer}: %.1fs | T_{Rnlds}: %.1fs | U: %.2f-%.2fm/s',dep_num,tau,tau_rd,round(vel_range(1),2),round(vel_range(2),2)),'fontsize',fs)
linkaxes(ax,'x')
xlim(ax(1),[0 10])

ylim(ax(1),[0 0.3]);
ylim(ax(2),[-2e-4 4.05e-4])

lgd = legend(ax(1),h20,{'ADV burst','spatial mean'},'fontsize',fs-2,'orientation','horizontal','location','northeast');

u_star = sqrt(mean(-uv_sgf,'omitnan'));
Cd = median(u_star^2/mean(vel_mag,'omitnan').^2);
fprintf('u*: %.4f m/s\nCd: %.4f\n',round(u_star,4),round(Cd,4))

%% save figures if desired
save_figs = 0;
if save_figs
    preamble = sprintf('dep%d/dep%d_%ds_%ds_%d-%d_',dep_num,dep_num,round(tau),round(tau_rd),round(100*vel_range(1)),round(100*vel_range(2)));
    print(figure(11),['stress_figures/' preamble 'velocity.png'],'-dpng','-r300')
    print(figure(12),['stress_figures/' preamble 'stress.png'],'-dpng','-r300')
    print(figure(13),['stress_figures/' preamble 'tke.png'],'-dpng','-r300')
end