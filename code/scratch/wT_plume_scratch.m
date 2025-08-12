% script to see if any w'T' signals show up in spider data
% 6 Mar 2023

try
    adcp;
catch
    load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\adcp.mat
    load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\rbr.mat
    load windows_singingflower_0824.mat
    load beam_distance_singingflower_0824.mat
end

nT = 0.01*[rbr.pos]/sqrt(2);
z = adcp.burst.range;
nlbls = cell(1,length(nT));
for i = 1:length(nT)
    nlbls{i} = sprintf('n = %.2fm',nT(i));
end

% window
j = 4;
tidx1 = adcp.burst.time>=windows(j,1) & adcp.burst.time<=windows(j,2);
tidx2 = rbr(1).time>=windows(j,1) & rbr(1).time<=windows(j,2);

% rbr
time = rbr(1).time(tidx2);
T = nan(length(nT),length(time));
for i = 1:length(nT)
    T(i,:) = rbr(i).values(tidx2);
end
figure(1); clf; hold on
plot(time,T)
legend(nlbls)
datetick
ylabel('T (C)')
grid
box
title(sprintf('spider temperature (section %d)',j))

% adcp
bn = 1;
adcp2 = struct();
adcp2.burst = struct();
adcp2.burst.time = adcp.burst.time(tidx1);
adcp2.burst.range = z;
for fld = {'vel','cor','amp'}
    adcp2.burst.(fld{1}) = adcp.burst.(fld{1})(tidx1,:,:);
end

% qa
adcp2.burst.vel(adcp2.burst.cor<0) = nan;
figure(2); clf
adcpQuickPlot(figure(2),adcp2,'vel',0.03*[-1 1],[0 Inf],[0 0.8],1);

% bins
z0 = beam_distance(j,bn);
zT = z0-nT;
dz = 0.03;

bin_ridx = nan(length(nT),2);
binlbls = cell(1,length(nT));
for i = 1:length(nT)
    bin_idx(i,:) = extrema(find(z>=(zT(i)-dz) & z<=(zT(i)+dz)));
    binlbls{i} = sprintf('n = %.2f +/- %.3fm',nT(i),dz);
end

% lagged cross correlations for velocities range-binned
bin_vel = nan(length(zT),length(adcp2.burst.time));
for i = 1:length(zT)
    bin_vel(i,:) = squeeze(mean(adcp2.burst.vel(:,bin_idx(i,1):bin_idx(i,2),bn),2,'omitnan'));
    bin_vel(i,:) = meanFilter(meanFilter(bin_vel(i,:),5),5);
end
figure(3); clf; hold on
plot(adcp2.burst.time,bin_vel);
datetick
ylabel('binned, mean-filtered beam vel (m/s)')
legend(binlbls)
grid
box
title(sprintf('beam %d, section %d',bn,j))

% plot time series first
figure(4); clf
clear ax
for i = 1:length(nT)
    ax(i) = axes(figure(4),'position',axgridpos(length(nT),1,i,0.1,0.06));
    plot(adcp2.burst.time,bin_vel(i,:)-mean(bin_vel(i,:),'omitnan'))
    ylim(0.03*[-1 1])
    ylabel('vel (m/s) [mean removed]')
    yyaxis right
    set(gca,'ycolor','k')
    plot(time,T(i,:)-mean(T(i,:),'omitnan'),'linewidth',1)
    ylim(1.5*[-1 1])
    ylabel('T (C) [mean removed]')
    datetick
    title(sprintf('%s (beam %d, section %d)',binlbls{i},bn,j))
    grid
    if i==1
        legend({'vel (blue)','T (orange)'},'location','best')
    end
end
linkaxes(ax,'x')

bin_vel(isnan(bin_vel)) = 0;

dt = 0.5;
k = ceil(200/dt);
tau = (-k:k)*dt;

rho_binned = nan(length(nT),2*k+1);
for i = 1:length(nT)
    rho_binned(i,:) = laggedCrossCorr(bin_vel(i,1:4:end),T(i,:),k);
end

figure(5); clf; hold on
plot(tau,rho_binned)
plot(tau,tau*0,'color',0.3*[1 1 1])
ylabel('cross-correlation (beam vel - temp)')
xlabel('lag (s)')
xlim(extrema(tau))
legend(binlbls)
title(sprintf('beam %d, section %d',bn,j))
ylim(max(abs(ylim))*[-1 1])
grid
box

%% save figs
fig_path = fullfile('../figures/weekly/20230306',sprintf('beam%d',bn));
save_figs = false;
if save_figs
    print(figure(1),fullfile(fig_path,'temp.png'),'-dpng','-r400');
    print(figure(3),fullfile(fig_path,'binned_smoothed_vel.png'),'-dpng','-r400');
    print(figure(4),fullfile(fig_path,'uT_timeseries.png'),'-dpng','-r400');
    print(figure(5),fullfile(fig_path,'uT_crossCorr.png'),'-dpng','-r400');
end

