% Script to make plots of ADV backscatter profiles over an entire
% deployment.
%
% KJW
% 8 Aug 2025

clear

raw_dir = 'F:/meltstake/data/raw';
fig_dir = 'F:/meltstake/figures/adv_backscatter';

% segment info
ms_tbl = loadMSInfo(23,'overview');
dep_names = ms_tbl.Folder;
ndeps = length(dep_names);

% plot setup
pad = [.03 .05 .1 .1];
shift = [0 0];
fs = 12;

for i = 1:ndeps
    % load adv data
    adv_file = fullfile(raw_dir,dep_names{i},'adv','adv.mat');
    if exist(adv_file)
        load(fullfile(raw_dir,dep_names{i},'adv','adv.mat'))
    else
        continue
    end
    fprintf('%02d: %s\n',i,dep_names{i})

    % set up figure and plot
    fig = figure(1); clf
    clear ax
    for j = 1:3
        ax(j) = axes(fig,'position',axgridpos(1,3,j,pad,shift),'fontsize',fs);
        pcolor(adv.pck_time,adv.pck_dist(2,:),adv.pck_amp(:,:,j)')
        shading flat
        ylim(ax(j),[15 350])
        text(0.01,1,sprintf('beam %d',j),'units','normalized','verticalalignment','top','fontsize',fs-1)
    end
    linkaxes(ax)

    % clim
    c_min = mean(adv.pck_amp(:,adv.pck_dist(2,:)>=30 & adv.pck_dist(2,:)<=40,:),'all');
    c_max = 0.9*max(extrema(adv.pck_amp));
    for j = 1:3
        clim(ax(j),[c_min c_max])
    end
    
    % colorbar
    cbar_pos = cbarpos(ax(3),.01,.02);
    cbar = colorbar(ax(3),'position',cbar_pos);

    % labels
    ylabel(ax(1),'distance [mm]','fontsize',fs)
    title(ax(2),strrep(dep_names{i},'_','\_'),'fontsize',fs)

    % save figure
    C = strsplit(dep_names{i},'_');
    print(fig,fullfile(fig_dir,sprintf('adv_backscatter_%s_%s_%s.png',C{2},C{3},C{1})),'-dpng','-r300')

end


