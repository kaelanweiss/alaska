% Script to perform melt rate calculations on meltstake deployments
%
% KJW
% 29 Jan 2024

clear

% load deployments
i = 28;
ms = loadADVPCK(i);
ms_tbl = loadMSInfo;

tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
seg_tbl = readtable(tbl_path,'sheet','manualwindows');
seg_tbl = seg_tbl(seg_tbl.Number==i,:);

load('F:/adv/mean_profile.mat')

%% trim a little bit
t0 = datetime(2024,7,12,20,20,40);
idxt = ms.pck_time >= t0;
ms.pck_time = ms.pck_time(idxt);
ms.pck_dist = ms.pck_dist(idxt,:);
ms.pck_amp = ms.pck_amp(idxt,:,:);

% mean profile
% ms.pck_amp = ms.pck_amp./repmat(mean_prof',[size(ms.pck_amp,1) 1 3]);

% try adding time-wise hanning filtering
for i = 1:size(ms.pck_amp,2)
    for j = 1:size(ms.pck_amp,3)
        ms.pck_amp(:,i,j) = hannFilter(ms.pck_amp(:,i,j),5);
    end
end

%% find ice edges
range_limit = [75 300];
[edges,tdrill] = findADVIceEdge(ms,ms.ndrill,ms.x0,ms.L,ms.L0,ms.hann_width,range_limit);

% plot
clear ax
figure(1000); clf
for j = 1:3
    ax(j) = axes(figure(1000),'position',axgridpos(1,3,j,[.025,.025,.08,.08]));
    hold on

    % pcolor
    pcolor(ms.pck_time,ms.pck_dist(2,:),ms.pck_amp(:,:,j)')
    shading flat
    clim(extrema(ms.pck_amp))
    
    % edge
    plot(edges.time,medianFilter(edges.pos(:,j),1),'k.')

end

ylabel(ax(1),'dist (mm)')
title(ax(1),strrep(ms.depname,'_','\_'))
cbar = colorbar(ax(3));
cbar.Position = cbarpos(gca,.003,.015);

linkaxes(ax)
xlim(extrema(ms.pck_time))
ylim(extrema(ms.pck_dist(2,:)))

%% save edges
proc_path = 'F:meltstake\data\proc';
rsp = input(sprintf('Save pck_edges.mat for %s [y/n]? ',ms.depname),'s');
if strcmpi(rsp,'y')
    fprintf('Saving %s\n',fullfile(proc_path,ms.depname,'adv','pck_edges.mat'))
    if ~exist(fullfile(proc_path,ms.depname,'adv'),'dir')
        mkdir(fullfile(proc_path,ms.depname,'adv'))
    end
    save(fullfile(proc_path,ms.depname,'adv','pck_edges.mat'),'edges','tdrill')
end

