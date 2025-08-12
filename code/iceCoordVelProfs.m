% Script to compute along-ice velocity profiles from Spider beams pointing
% in vertical directions.
%
% KJW
% 17 Feb 2023

%clearvars -except fig

% load data
try
    adcp;
catch
    load F:/alaska2022/data/iceberg_surveys/mat/20220824_singingflower/spider/adcp.mat
    load windows_singingflower_0824.mat
    load beam_distance_singingflower_0824.mat
end

% beam angle
theta = 25;

% coordinate transform from instrument to ice-normal coordinates
n_trans = @(z,zw,theta,beta) (-z+zw*(tand(theta)*sind(beta)*cosd(beta)+1))/cosd(beta);

% calculate ice slope (beta)
zw2 = beam_distance(:,2);
zw4 = beam_distance(:,4);
beta = atand((zw2-zw4)*cotd(theta)./(zw2+zw4));

% calculate beam/along-ice separation
phi2 = 90-theta-beta;
phi4 = 90-theta+beta;

% choose window
j = 4;
idxt = adcp.burst.time>=windows(j,1) & adcp.burst.time<=windows(j,2);

% qa a little bit
qa_cor = adcp.burst.cor(idxt,:,:)<50;
vel = adcp.burst.vel(idxt,:,:);
vel(qa_cor) = nan;
adcp.burst.vel(idxt,:,:) = vel;
vel_mean = squeeze(mean(adcp.burst.vel(idxt,:,:),'omitnan'));
b2 = vel_mean(:,2);
b4 = vel_mean(:,4);

% plot panels
rmax = 0.9;
% adcpQuickPlot(fig(1),adcp,'vel',0.04*[-1 1],windows(j,:),[0 rmax],1);
% adcpQuickPlot(fig(2),adcp,'cor',[50 100],windows(j,:),[0 rmax],1);
% adcpQuickPlot(fig(3),adcp,'amp',[65 88],windows(j,:),[0 rmax],1);
% figure(1)

% transform coordinates
n2 = n_trans(adcp.burst.range,zw2(j),theta,beta(j));
n4 = n_trans(adcp.burst.range,zw4(j),theta,beta(j));

nmin = -0.2;
idxn2 = n2>=nmin;
idxn4 = n4>=nmin;

% transform velocities
w2 = b2/cosd(phi2(j));
w4 = -b4/cosd(phi4(j));

%% plot
lw = 0.9;
figure(4); clf
axes(figure(4),'position',axgridpos(1,2,1,.08,.12))
hold on
plot(n2(idxn2),b2(idxn2),'.-','linewidth',lw)
plot(n4(idxn4),b4(idxn4),'.-','linewidth',lw)
xlabel('n (m)')
ylabel('beam vel (m/s)')
grid
box on
legend({'b2 (up)','b4 (dn)'},'location','best')
title(sprintf('SF0824 Section%d \\beta=%.1f^\\circ',j,beta(j)))

axes(figure(4),'position',axgridpos(1,2,2,.08,.12))
hold on
plot(n2(idxn2),w2(idxn2),'.-','linewidth',lw)
plot(n4(idxn4),w4(idxn4),'.-','linewidth',lw)
xlabel('n (m)')
ylabel('scaled along-ice vel (m/s)')
grid
box on
legend({'b2 (up)','b4 (dn)'})

%print(figure(4),fullfile('../figures/weekly/20230220',sprintf('ice_coord_vert_vel_section%d.png',j)),'-dpng','-r400');