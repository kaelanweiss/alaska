% Script to plot processed data from each meltstake segment for easy
% inspection.
%
% KJW
% 7 Apr 2025
clear

proc_dir = 'F:/meltstake/data/proc';
fig_dir = 'F:/meltstake/figures/segment_summary';

% use bin-mapped data
use_binmap = 0;

% segment info
ms_tbl = loadMSInfo('segments');

dep_nums = unique(ms_tbl.Number);
ndeps = length(dep_nums);

% plotting
lw = 1;
fs = 12;
figsize = [15 18];
pad = [.05 .02 .15 .07];
shift = [-.05 0];
vel_comp_order = [1 3 2];
vel_lbls = {'u','v (x2)','w','U_{max}'};
T_clr = [0.2 1 0.5]/1.2;
T_clrs = T_clr'*[1 0.5 0];
T_lbls = {'near','mid','far'};

% loop through deployment
for i = 19:ndeps
    nsegs = sum(ms_tbl.Number==dep_nums(i));
    dep_name = ms_tbl.Folder{ms_tbl.Number==dep_nums(i) & ms_tbl.Window==1};
    for j = 1:nsegs
        % load data
        if use_binmap
            load(fullfile(proc_dir,dep_name,'adcp_bin_map',sprintf('adcp%d.mat',j)))
        else
            load(fullfile(proc_dir,dep_name,sprintf('adcp%d.mat',j)))
        end
        load(fullfile(proc_dir,dep_name,sprintf('T%d.mat',j)))

        % mean velocity components
        adcp.burst.vel_ice(:,:,3) = 2*adcp.burst.vel_ice(:,:,3);
        vel_mean = squeeze(mean(adcp.burst.vel_ice(:,:,vel_comp_order),2,'omitnan'));
        for k = 1:3
            vel_mean(:,k) = hannFilter(vel_mean(:,k),2*adcp.burst.samplerate);
        end
        u_max = hannFilter(adcp.burst.u_max,2*adcp.burst.samplerate);

        % figure
        fig = figure(1); clf
        setFigureSize(fig,figsize);
        clear ax
        
        % velocity pcolors
        vel_lim = min([max(abs(extrema(adcp.burst.vel_ice))) 4*std(adcp.burst.vel_ice(:,:,vel_comp_order),0,'all','omitnan')]);
        for k = 1:3
            ax(k) = axes(fig,'position',axgridpos(5,1,k,pad,shift));
            pcolor(ax(k),adcp.burst.time,adcp.burst.range,adcp.burst.vel_ice(:,:,vel_comp_order(k))')
            shading flat
            cmocean('bal')
            clim(vel_lim*[-1 1])
            text(0.02,0.93,vel_lbls{k},'fontsize',fs,'units','normalized')
            ax(k).FontSize = fs-2;
        end
        cbar = colorbar(ax(3),'position',cbarpos(ax(3),0.015,.03));
        cbar.Position(4) = sum(ax(1).Position([2 4])) - ax(3).Position(2);


        % velocity components
        clear hvel
        ax(4) = axes(fig,'position',axgridpos(5,1,4,pad,shift),'fontsize',fs-2);
        hold on
        box on
        plot(adcp.burst.time([1 end]),[0 0],'color',0.4*[1 1 1])
        for k = 1:3
            hvel(k) = plot(adcp.burst.time,vel_mean(:,k),'linewidth',lw,'color',colors(k));
        end
        hvel(4) = plot(adcp.burst.time,u_max,'k','linewidth',lw);
        legend(hvel,vel_lbls,'location','eastoutside','fontsize',fs-2)
        ax(4).Position = axgridpos(5,1,4,pad,shift);

        % temperature
        ax(5) = axes(fig,'position',axgridpos(5,1,5,pad,shift),'fontsize',fs-2);
        hold on
        box on
        for k = 1:size(T.T,2)
            plot(T.time,T.T(:,k),'linewidth',lw,'color',T_clrs(:,k));
        end
        legend(T_lbls,'location','eastoutside','fontsize',fs-2)
        ax(5).Position = axgridpos(5,1,5,pad,shift);
    
        % clean up axes
        for k = 1:4
            ax(k).XTickLabels = [];
        end
        linkaxes(ax,'x')
        xlim(ax(1),adcp.burst.time([1 end]))
        % labels
        for k = 1:3
            ylabel(ax(k),'range [m]','fontsize',fs)
        end
        ylabel(ax(4),'mean vel [m/s]','fontsize',fs)
        ylabel(ax(5),'T_w [\circC]','fontsize',fs)
        cbar.Label.String = 'velocity [m/s]';
        cbar.Label.FontSize = fs;

        % title
        depth = ms_tbl.depth(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        dur = ms_tbl.Duration(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        m_obs = msTable2Vector(ms_tbl.m(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j))*.24;
        if use_binmap
            title(ax(1),sprintf('dep %d seg %d | %.1f m, %d min, %.1f m/day (binmap)',dep_nums(i),j,round(depth,1),round(dur),round(m_obs,1)),'fontsize',fs)
        else
            title(ax(1),sprintf('dep %d seg %d | %.1f m, %d min, %.1f m/day',dep_nums(i),j,round(depth,1),round(dur),round(m_obs,1)),'fontsize',fs)
        end

        % save figure
        if use_binmap
            print(fig,fullfile(fig_dir,'adcp_bin_map',sprintf('segment_summary_%02d_%02d.png',dep_nums(i),j)),'-dpng','-r300')
        else
            print(fig,fullfile(fig_dir,sprintf('segment_summary_%02d_%02d.png',dep_nums(i),j)),'-dpng','-r300')
            a = 1;
        end
    end
end
