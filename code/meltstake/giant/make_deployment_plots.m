% script to make quick plots of the meltstake deployment(s)
%
% KJW
% 9 Aug 2026

clear

depname = '20260730_1530_ms01';
raw_dir = 'F:\GIANT2026\meltstake\data\raw';
proc_dir = 'F:\GIANT2026\meltstake\data\proc';
fig_dir = 'F:\GIANT2026\meltstake\figures';

fs = 11;
lw = 1;

%% load T and S things
% T
load(fullfile(raw_dir,depname,'solo','T.mat'))

% MS CTD
load(fullfile(raw_dir,depname,'ctd_ms','ctd_ms.mat'))

% ROV CTD
load(fullfile(raw_dir,depname,'ctd_rov','ctd_rov.mat'))


%% plot small files
% T
figure(1); clf
hold on
for i = 1:length(T)
    plot(T(i).time,T(i).values)
end
plot(ctd_ms.time,ctd_ms.T,'k-')
box on
ylabel('temperature [^\circC]','fontsize',fs)
legend({'10cm','30cm','50cm','CTD (30cm)'},'fontsize',fs-1,'orientation','horizontal')

% S
figure(2); clf
plot(ctd_ms.time,ctd_ms.S,'k-','linewidth',lw)
box on
ylabel('salinity [psu]','fontsize',fs)

%% ROV CTD
figure(3); clf
for i = 1:2
    ax3(i) = subplot(1,2,i);
    hold on
    axis ij
    box on
end
ylabel(ax3(1),'depth [m]','fontsize',fs)
xlabel(ax3(1),'temperature [^\circC]','fontsize',fs)
xlabel(ax3(2),'salinity [psu]','fontsize',fs)

% deploy cast
t1 = datetime(2026,8,3,14,15,20);
t2 = datetime(2026,8,3,14,36,53);
idxt = ctd_rov.time > t1 & ctd_rov.time < t2;
plot(ax3(1),ctd_rov.values(idxt,2),ctd_rov.values(idxt,4),'linewidth',lw)
plot(ax3(2),ctd_rov.values(idxt,5),ctd_rov.values(idxt,4),'linewidth',lw)
% outer cast
t1 = datetime(2026,8,3,14,57,45);
t2 = datetime(2026,8,3,15,4,19);
idxt = ctd_rov.time > t1 & ctd_rov.time < t2;
plot(ax3(1),ctd_rov.values(idxt,2),ctd_rov.values(idxt,4),'linewidth',lw)
plot(ax3(2),ctd_rov.values(idxt,5),ctd_rov.values(idxt,4),'linewidth',lw)

legend(ax3(1),{'deploy cast','outer cast'},'fontsize',fs-1,'location','southeast')

%% adcp
load(fullfile(raw_dir,depname,'adcp','adcp.mat'))
adcp.burst.vel_uvw = adcp.burst.vel_ice(:,:,[1 3 2]);
%% adcp plots
ax4 = adcpQuickPlot(figure(4),adcp,'vel_uvw',.2,NaT,[0 2],10);
plbls4 = addPanelLabels(ax4,fs+1);
plbls4(1).String = 'u (along-ice)';
plbls4(2).String = 'v (ice-normal)';
plbls4(3).String = 'w (vertical)';

%% adcp average
idxr = adcp.burst.range < 0.4;
uvw = squeeze(mean(adcp.burst.vel_uvw(:,idxr,:),2,'omitnan'));
uw_mag = vecnorm(uvw(:,[1 3]),2,2);

% smooth 
uvw_smth = uvw;
for i = 1:3
    uvw_smth(:,i) = hannFilter(uvw_smth(:,i),8*30);
end
uw_mag_smth = hannFilter(uw_mag,8*30);

figure(5); clf
hold on
plot(adcp.burst.time,uw_mag_smth,'k-','linewidth',2*lw)
for i = 1:3
    plot(adcp.burst.time,uvw_smth(:,i),'linewidth',lw,'color',colors(i))
end

box on
ylabel('velocity [m/s]','FontSize',fs)
legend({'sqrt(u^2+w^2)','u','v','w'},'fontsize',fs-1,'orientation','horizontal')
xlim(extrema(adcp.burst.time))

%% psd
uvw(isnan(uvw)) = 0;
[P,f] = psd_matlab(uvw,adcp.burst.samplerate,'notaper');

% binning setup
nf = [4 8 16 32];
nf_bins = [0 1e-2 1e-1 1e0];
[Pbin,fbin] = progressiveBin(f,P(:,1),nf,nf_bins);
for i = 2:3
    Pbin(:,i) = progressiveBin(f,P(:,i),nf,nf_bins);
end

figure(6); clf
plot(fbin,Pbin,'linewidth',lw)
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('f [cps]','fontsize',fs)
ylabel('PSD [m^2 s^{-2} cps^{-1}]','fontsize',fs)
legend({'u','v','w'},'fontsize',fs-1)

%% adv
load(fullfile(raw_dir,depname,'adv','adv.mat'))
load(fullfile(proc_dir,depname,'adv','pck_edges.mat'))

%% adv distance
figure(7); clf
clear ax7
dt = 10;
for i = 1:3
    ax7(i) = axes(figure(7),'position',axgridpos(1,3,i,[0.02,0.05,0.1,0.12]));
    box on
    pcolor(adv.pck_time(1:dt:end),adv.pck_dist(1,:),adv.pck_amp(1:dt:end,:,i)')
    shading flat
    clim(extrema(adv.pck_amp))
    hold on
    % plot(edges.time,edges.pos(:,i)-10,'k-')
    if i == 1
        ylabel('dist (mm)')
    end
    clim([50 200])
end
cbar = colorbar;
cbar.Position = [0.65+.25+.01 .12 .02 .76];

linkaxes(ax7)
ylim(ax7(1),[0 450])