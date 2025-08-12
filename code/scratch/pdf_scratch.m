% script to play around with pdfs of spider velocity
% 3 Mar 2023

% load
try
    adcp
catch
    load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\adcp.mat
    load windows_singingflower_0824.mat
end

% window
j = 4;
idxt = adcp.burst.time>=windows(j,1) & adcp.burst.time<=windows(j,2);

vel = adcp.burst.vel(idxt,:,:);
cor = adcp.burst.cor(idxt,:,:);
amp = adcp.burst.amp(idxt,:,:);

% qc
vel_qc = vel;
vel_qc(cor<0) = nan;

%% histogram
beam = 3;
rmax = 0.68;
[~,J] = min(abs(adcp.burst.range-rmax));

bmax = 0.05;
bstep = 0.0025;
bins = -bmax:bstep:bmax;

x = bins(1:end-1) + 0.5*diff(bins);

N = nan(J,length(bins)-1);
for i = 1:J
    N(i,:) = histcounts(vel_qc(:,i,beam),bins)/length(find(~isnan(vel_qc(:,i,beam))));
end

% plot
scale = 0.75;
figure(10+beam); clf; hold on
for i = J:-1:1
    area(x,N(i,:)*scale+adcp.burst.range(i),adcp.burst.range(i),'facecolor',0.6*[1 1 1])
    area(x,0*x+adcp.burst.range(i),'FaceColor','w')
    plot(x,0*x+adcp.burst.range(i),'k:')
    plot(x,N(i,:)*scale+adcp.burst.range(i),'k-')
end

zline = [adcp.burst.range(1) adcp.burst.range(J)+1.2*scale*max(N(J,:))];
plot([0 0],zline,'-','color',0.3*[1 1 1])
title(sprintf('beam%d',beam))

yticks = adcp.burst.range(1:2:J);
set(gca,'ytick',yticks)
set(gca,'yticklabel',strsplit(sprintf('%.2f,',yticks),','))
ylim(zline)