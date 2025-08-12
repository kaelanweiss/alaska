% Script to plot CTD profiles for meltstake deployments. Profiles need to
% be collected into ctd.mat files in the processed deployment folders.
%
% NOTE: currently set up to work on the 2024 deployment ctd file format,
% deployments from 2023 will break the code because of different field
% names
%
% KJW
% 5 March 2024

clear

% load deployment data
ms_tbl = loadMSInfo;

% directories
fig_dir = 'F:meltstake/figures';
proc_dir = 'F:meltstake\data\proc';

% choose deployments
deps = 20:28;

% plot formatting
pad = [.08 .05 .12 .12];
shift = [0 0];
lw = 1.5;
fs = 11;
clrs = [1 7];

% save figures?
save_figs = 0;

%% loop through deployments
for i = deps
    % load ctd data
    load(fullfile(proc_dir,ms_tbl.Folder{i},'ctd.mat'))
    
    % deploment depth
    ms_z = str2double(ms_tbl.Depth_m_{i});

    % figure and axes
    fig = figure(i); clf
    setFigureSize(fig,[16 10]);
    clear ax
    for j = 2:-1:1
        ax(j) = axes(fig,'position',axgridpos(1,2,j,pad,shift),'fontsize',fs-2);
        hold on
        axis ij
        box on
        plot(ax(j),[0 40],ms_z*[1 1],'k--','linewidth',1)
    end
    
    % plot
    % T
    TLIM = nan(length(ctd),2);
    for j = 1:length(ctd)
        if ~isempty(ctd(j).T)
            plot(ax(1),ctd(j).T,ctd(j).depth,'color',0.5*[1 1 1])
            plot(ax(1),mean(ctd(j).T,2,'omitnan'),ctd(j).depth,'color',colors(clrs(j)),'linewidth',lw)
            TLIM(j,:) = extrema(ctd(j).T);
        end
    end

    % S
    clear h
    SLIM = nan(length(ctd),2);
    for j = 1:length(ctd)
        if ~isempty(ctd(j).T)
            plot(ax(2),ctd(j).S,ctd(j).depth,'color',0.5*[1 1 1])
            plot(ax(2),mean(ctd(j).S,2,'omitnan'),ctd(j).depth,'color',colors(clrs(j)),'linewidth',lw)
            SLIM(j,:) = extrema(ctd(j).S);
        end
        h(j) = plot(ax(2),0,0,'color',colors(clrs(j)),'linewidth',lw);
    end

    % legend
    lgd_lbls = cell(1,length(ctd));
    for j = 1:length(ctd)
        lgd_lbls{j} = sprintf('%s (%d)',ctd(j).vessel,length(ctd(j).time));
    end
    legend(ax(2),h,lgd_lbls,'location','southwest','fontsize',fs-1)

    % labels
    xlabel(ax(1),'temperature [\circC]','fontsize',fs)
    xlabel(ax(2),'salinity [psu]','fontsize',fs)
    ylabel(ax(1),'depth [m]','fontsize',fs)

    % titles
    t1 = ms_tbl.Start(i);
    t2 = ms_tbl.End(i);
    t1.Format = 'uuuu-MM-dd HH:mm';
    t2.Format = 'HH:mm';

    title(ax(1),sprintf('%s - %s UTC',t1,t2),'fontsize',fs)
    title(ax(2),sprintf('meltstake depth: %d m',round(ms_z)),'fontsize',fs)

    % limits
    TLIM = [min(floor(TLIM(:,1))) max(ceil(TLIM(:,2)))];
    SLIM = [min(floor(SLIM(:,1))) max(ceil(SLIM(:,2)))];
    SLIM(1) = max([SLIM(1) 20]);
    xlim(ax(1),TLIM)
    xlim(ax(2),SLIM)

    % save figure
    if save_figs
        if ~exist(fullfile(fig_dir,'fjord_ctd'),'dir')
            mkdir(fullfile(fig_dir,'fjord_ctd'))
        end
        fname_profs = sprintf('ctd_profiles_%s.png',ms_tbl.Folder{i});
        print(fig,fullfile(fig_dir,'fjord_ctd',fname_profs),'-dpng','-r350')
    end

end

