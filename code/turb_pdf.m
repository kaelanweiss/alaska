% script to calculate flow pdfs of spider velocity
% 8 Mar 2023

% load
try
    adcp
catch
    load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\adcp.mat
    load windows_singingflower_0824.mat
end

[~,nc,nb] = size(adcp.burst.vel);
nw = size(windows,1);

%% histogram
bmax = 0.08;
bstep = 0.0025;
bins = -bmax:bstep:bmax;

x = bins(1:end-1) + 0.5*diff(bins);
nx = length(x);
pdf = nan(nx,nc,nb);

clear vel_counts
vel_counts(nw) = struct();

for k = 1:nw
    % window
    idxt = adcp.burst.time>=windows(k,1) & adcp.burst.time<=windows(k,2);

    vel_counts(k).x = x;
    vel_counts(k).pdf = pdf;

    vel = adcp.burst.vel(idxt,:,:);
    cor = adcp.burst.cor(idxt,:,:);

    % qc
    vel_qc = vel;
    vel_qc(cor<50) = nan;

    for j = 1:nb
        % beam
        for i = 1:nc
            % cell
            vel_counts(k).pdf(:,i,j) = histcounts(vel_qc(:,i,j),bins)/length(find(~isnan(vel_qc(:,i,j))));
        end
    end
end

%% plot
k = 6; % window
j = 4; % beam
z0 = beam_distance(k,j);

pdf = vel_counts(k).pdf;
imax = floor((z0-0.12)/0.02)-2;
scale = 0.75;
figure(10); clf; hold on
for i = imax:-1:1
    area(100*x,pdf(:,i,j)*scale+adcp.burst.range(i),adcp.burst.range(i),'facecolor',0.3*[1 1 1])
    area(100*x,0*x+adcp.burst.range(i),'FaceColor','w')
    plot(100*x,0*x+adcp.burst.range(i),'k:')
    plot(100*x,pdf(:,i,j)*scale+adcp.burst.range(i),'k-','linewidth',0.9)
end

zline = [adcp.burst.range(1) adcp.burst.range(imax)+1.2*scale*max(pdf(:,imax,j))];
plot([0 0],zline,'-','color',0.5*[1 1 1])
title(sprintf('velocity histogram (beam%d, section%d)',j,k))

yticks = adcp.burst.range(flip(imax:-2:1));
set(gca,'ytick',yticks)
set(gca,'yticklabel',strsplit(sprintf('%.0f,',-100*(yticks-z0)),','))
set(gca,'xtick',-100*bmax:2:100*bmax)
ylim(zline)
box
xlabel('beam velocity (cm/s)')
ylabel('distance from ice (cm)')

if true
    print(figure(10),sprintf('../figures/corvallis-march23/pdfs/s%d_b%d_SF0824.png',k,j),'-dpng','-r600')
end