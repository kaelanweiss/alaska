%clear

load F:/meltstake/data/raw/ms02_20240613_2230/adv/adv.mat

%%
amp = adv.pck_amp;
%amp(repmat(adv.pck_dist,[1 1 3])<40) = nan;
pad1 = [0.02,0.05,0.1,0.12];
figure(11); clf
clear ax
for i = 1:3
    ax(i) = axes(figure(11),'position',axgridpos(1,3,i,pad1));
    box on
    pcolor(repmat(adv.pck_time,[1 size(adv.pck_dist,2)]),adv.pck_dist,amp(:,:,i))
    %pcolor(repmat(adv.pck_time,[1 275]),repmat(1:275,[length(adv.pck_time) 1]),amp(:,:,i))
    shading interp
    clim(extrema(adv.pck_amp))
    hold on
%     plot(adv.hdr_time,(medianFilter(adv.dp(:,i,1),1)+50)/sind(60),'k-')
%     plot(adv.hdr_time,(medianFilter(adv.dp(:,i,1),1)+50)/sind(60),'k-')
%     plot(adv.hdr_time,adv.dp(:,i,1),'k-')
%     plot(adv.hdr_time,adv.dp(:,i,2),'k--')
%     plot(adv.hdr_time,adv.dprobe_ave,'r')
%     plot(adv.hdr_time,adv.dsvol_ave(:,1)+150,'b')
    if i == 1
        ylabel('dist (mm)')
    end
end
cbar = colorbar;
cbar.Position = [0.65+.25+.01 .12 .02 .76];

linkaxes(ax)

%%
pad2 = [0.08,0.05,0.1,0.12];
figure(12); clf
axes(figure(12),'position',axgridpos(1,2,1,pad2))
box on
pcolor(repmat(adv.pck_time,[1 275]),adv.pck_dist,amp(:,:,i))
%pcolor(repmat(adv.pck_time,[1 275]),repmat(1:275,[length(adv.pck_time) 1]),amp(:,:,i))
shading flat
clim(extrema(adv.pck_amp))
ylabel('dist (mm)')

axes(figure(12),'position',axgridpos(1,2,2,pad2))
box on
%pcolor(repmat(adv.pck_time,[1 275]),adv.pck_dist,amp(:,:,i))
pcolor(repmat(adv.pck_time,[1 275]),repmat(1:275,[length(adv.pck_time) 1]),amp(:,:,i))
%pcolor(adv.pck_time,adv.pck_dist(1,:),amp(:,:,i)')
shading flat
clim(extrema(adv.pck_amp))
ylabel('cell number')


cbar = colorbar;
cbar.Position = [0.65+.25+.01 .12 .02 .76];

%%
amp_all = prod(adv.pck_amp/max(adv.pck_amp,[],'all'),3);

% axes(figure(1),'position',axgridpos(1,4,4,0.02,0.05,0.05,0.08))
% box on
% pcolor(repmat(adv.pck_time,[1 275]),adv.pck_dist,amp_all)
% shading flat

figure(2); clf
surf(adv.pck_time,adv.pck_dist(1,:),adv.pck_amp(:,:,1)')
shading flat
ylabel('dist (mm)')

amp_smooth = adv.pck_amp;
for i = 1:size(amp_smooth,2)
    for j = 1:size(amp_smooth,3)
        amp_smooth(:,i,j) = hannFilter(amp_smooth(:,i,j),7);
    end
end

figure(3); clf
surf(adv.pck_time,adv.pck_dist(1,:),amp_smooth(:,:,1)')
shading flat
ylabel('dist (mm)')