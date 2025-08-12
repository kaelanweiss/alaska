
try
    adcp;
catch
    load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\adcp.mat
    load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\rbr.mat
    load windows_singingflower_0824.mat
    load beam_distance_singingflower_0824.mat
end

j = 6;
imax = 40;

beams = [3 4];

idxt = adcp.burst.time>=windows(j,1) & adcp.burst.time<=windows(j,2);

z = repmat(adcp.burst.range(1:imax)',[1 length(beams)]);
for i = 1:length(beams)
    z(:,i) = -100*(z(:,i)-beam_distance(j,beams(i)))*mean(beam_distance(j,beams))/beam_distance(j,beams(i))+1;
end

vel = adcp.burst.vel(idxt,1:imax,beams);
cor = adcp.burst.cor(idxt,1:imax,beams);

vel(cor<50) = nan;

vel_mean = 100*squeeze(mean(vel,'omitnan'));

for i = 1:length(beams)
    [maxv,imaxv] = max(abs(vel_mean(:,i)));
    vel_mean(:,i) = (vel_mean(imaxv,i)/maxv)*(1/sind(25))*vel_mean(:,i);
end

rbr_dist = nan(1,length(rbr));
T_mean = nan(1,length(rbr));
T_std = nan(1,length(rbr));
for i = 1:length(rbr)
    idxt_rbr = rbr(i).time>=windows(j,1) & rbr(i).time<=windows(j,2);
    rbr_dist(i) = rbr(i).pos/sqrt(2);
    T_mean(i) = mean(rbr(i).values(idxt_rbr));
    T_std(i) = std(rbr(i).values(idxt_rbr));
end

%% plot
rmin = 4;
ls = {'.-','.:','.:'};
ms = 8;
lw = 1;
figure(40+j); clf; hold on
%set(figure(40+j),'position',[-1500 130 520 331])
for i = 1:length(beams)
    idxr = z(:,i) > rmin;
    plot(z(idxr,i),vel_mean(idxr,i),ls{i},'linewidth',lw,'color','k','markersize',ms)
end
grid
box
xlim([0 55])
xlabel('distance from ice (cm)')
ylabel('vel (cm/s)')
ylim([-1 5])

yyaxis right
set(gca,'YColor','r')
%plot([0 rbr_dist(1)],[0 T_mean(1)],'r.--','linewidth',0.8)
errorbar(rbr_dist,T_mean,T_std,'r.-','linewidth',lw,'markersize',ms)
ylabel('temp (C)')
ylim([4 7])

legend({'horizontal flow','vertical flow'},'location','northwest')
