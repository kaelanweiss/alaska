% script to look at the effect of mean flow on BL thermal structure

try
    adcp;
catch
    load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\adcp.mat
    load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\rbr.mat
end

pos = round([rbr.pos]/sqrt(2));
poslbls = cell(size(pos));
for i = 1:length(pos)
    poslbls{i} = sprintf('n=%dcm',pos(i));
end

% filtering window
k = ceil(8*60*5);

% create copy of velocity for filtering
vel2 = adcp.burst.vel;
vel2(adcp.burst.cor<70) = nan;

% get rid of data where ADCP was moving (judged by in-ice velocity)
vel_ice = squeeze(mean(adcp.burst.vel(:,35:55,:),2));
for i = 1:4
    [~,idxi] = sigmaFilter(vel_ice(:,i),1,1,1);
    vel_temp = vel2(:,:,i);
    vel_temp(repmat(idxi,[1 size(vel_temp,2)])) = nan;
    vel2(:,:,i) = vel_temp;
end

% collect temperature data
T = nan(length(rbr(1).time),length(rbr));
Tmean = nan*T;
for i = 1:length(rbr)
    T(:,i) = rbr(i).values;
    Tmean(:,i) = hannFilter(T(:,i),round(k/4)+1);
end

% run filter
tic;
for j = 1:4
    for i = 1:85
        %vel2(:,i,j) = meanFilter(meanFilter(vel2(:,i,j),k,'omitnan'),k,'omitnan');
        vel2(:,i,j) = hannFilter(vel2(:,i,j),k);
        fprintf('%d: %2d\n',j,i)
    end
end
t_elapsed = toc;

% mean outer velocity
vel_outer = squeeze(mean(vel2(:,1:20,[1 3]),2,'omitnan'));
vel_outer(:,2) = -vel_outer(:,2);

% plot velocity
adcp2 = adcp;
adcp2.burst.vel = vel2;
adcpQuickPlot(figure(1),adcp2,'vel',0.03*[-1 1],[0 Inf],[0 0.9],1);

% plot temp
figure(2); clf; hold on
plot(rbr(1).time,Tmean,'linewidth',0.9)
xlim(extrema(rbr(1).time))
datetick('x','mmmdd HH:MM','keeplimits')
ylabel('temperature (C)')
legend(poslbls,'location','best')
grid
box

% plot mean outer velocity
figure(3); clf; hold on
bmlbls = {'beam1','-beam3'};
plot(adcp.burst.time,vel_outer,'linewidth',0.9)
xlim(extrema(adcp.burst.time))
datetick('x','keeplimits')
ylabel({'beam vel (m/s)','averaged [0.12,0.5]m'})
legend(bmlbls,'location','eastoutside')
grid
box
