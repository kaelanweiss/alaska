% script to plot summaries of a few segments
%
% KJW
% 19 Apr 2025

clear

addpath('..')

proc_dir = 'F:/meltstake/data/proc';
load clrs_egu.mat

% segment info
ms_tbl = loadMSInfo(26:28,'manualwindows');

dep_nums = unique(ms_tbl.Number);
ndeps = length(dep_nums);

% choose segment(s) from each deployment
seg_nums = {[2], [5], [6]};

% plotting
lw = 1.2;
fs = 13.5;
figsize = [16 18];
pad = [.05 .025 .18 .07];
shift = [-.05 0];
vel_comp_order = [1 3 2];
vel_lbls = {'u','v','w','U_{max}'};
T_clr = [0.2 1 0.5]/1.2;
T_clrs = T_clr'*[1 0.5 0];
T_lbls = {'0.1 m','0.3 m','0.5 m'};
prof_ylims = [-.05 .03;-.05 .03;-.01 .06];

% loop through deployment
for i = 1:ndeps
    dep_name = ms_tbl.Folder{ms_tbl.Number==dep_nums(i) & ms_tbl.Window==1};
    for j = seg_nums{i}
        % load data
        load(fullfile(proc_dir,dep_name,sprintf('adcp%d.mat',j)))
        load(fullfile(proc_dir,dep_name,sprintf('T%d.mat',j)))

        % mean velocity components
        vel_mean = squeeze(mean(adcp.burst.vel_ice(:,:,vel_comp_order),2,'omitnan'));
        for k = 1:3
            vel_mean(:,k) = hannFilter(vel_mean(:,k),2*adcp.burst.samplerate);
        end
        u_max = hannFilter(adcp.burst.u_max,2*adcp.burst.samplerate);

        % velocity magnitude
        vel_mag = mean(vecnorm(adcp.burst.vel_ice,2,3),'omitnan')';
        std_mag = std(vecnorm(adcp.burst.vel_ice,2,3),'omitnan')';

        % figure
        fig = figure(i); clf
        setFigureSize(fig,figsize);
        clear ax
        
        % velocity pcolors
        vel_lim = min([max(abs(extrema(adcp.burst.vel_ice))) 4*std(adcp.burst.vel_ice(:,:,vel_comp_order),0,'all','omitnan')]);
        rmax = ms_tbl.rmax(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        y = rmax/cosd(25) - adcp.burst.range;
        for k = 1:3
            ax(k) = axes(fig,'position',axgridpos(5,1,k,pad,shift));
            pcolor(ax(k),adcp.burst.time,y,adcp.burst.vel_ice(:,:,vel_comp_order(k))')
            shading flat
            cmocean('bal')
            clim(vel_lim*[-1 1])
            text(0.02,0.93,vel_lbls{k},'fontsize',fs,'units','normalized','backgroundcolor','none')
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
            hvel(k) = plot(adcp.burst.time,vel_mean(:,k),'linewidth',lw,'color',vel_clrs(k,:));
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
        ylim(ax(4),0.4*[-1 1])
        ylim(ax(5),[3.5 7])
        set(ax(4),'clipping','off')
        % labels
        for k = 1:3
            ylabel(ax(k),'y [m]','fontsize',fs)
            ylim(ax(k),[0 0.515])%max(y)])
            clim(ax(k),0.35*[-1 1])
            set(ax(k),'ytick',0:0.2:1)
        end
        ylabel(ax(4),'mean vel [m/s]','fontsize',fs)
        ylabel(ax(5),'T_w [\circC]','fontsize',fs)
        cbar.Label.String = 'velocity [m/s]';
        cbar.Label.FontSize = fs;

        % title
        depth = ms_tbl.depth(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        dur = ms_tbl.Duration(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        m_obs = msTable2Vector(ms_tbl.m(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j))*.24;
%         title(ax(1),sprintf('dep %d seg %d | %.1f m, %d min, m_{obs}=%.1f m/day',dep_nums(i)-25,j,round(depth,1),round(dur),round(m_obs,1)),'fontsize',fs)
        
        % velocity profiles
        fig_prof = figure(i+3); clf
        setFigureSize(fig_prof,[11 6])
        axp = axes(fig_prof,'position',axgridpos(1,1,1,[0 0 .1 .15],[0.055 0.06]),'fontsize',fs-2);
%         axp = axes(fig_prof,'fontsize',fs-2);
        hold on
        clear hprof
        % compute profile things
        vel_prof = squeeze(mean(adcp.burst.vel_ice(:,:,vel_comp_order),1,'omitnan')); % component means
        mag_comps = vecnorm(vel_prof,2,2); % magnitude of mean components
        tau = ms_tbl.tau_decorr(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j); % decorrelation
        Neff = 60*ms_tbl.Duration(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j)/ms_tbl.tau_decorr(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        prof_std = squeeze(std(adcp.burst.vel_ice(:,:,vel_comp_order),0,1,'omitnan'))/sqrt(Neff); % component standard error
        dmag_comps = sqrt(sum((vel_prof.*prof_std).^2,2))./mag_comps;
        
        plot(axp,[0 1],[0 0],'--','color',0.*[1 1 1],'linewidth',lw)
        for k = 1:3
            plotProfConf(axp,y',vel_prof(:,k)-prof_std(:,k),vel_prof(:,k)+prof_std(:,k),vel_clrs(k,:),0.3);
        end
        plotProfConf(axp,y',mag_comps-dmag_comps,mag_comps+dmag_comps,0.5*vel_clrs(k,:),0.3);
        plotProfConf(axp,y',vel_mag-std_mag/sqrt(Neff),vel_mag+std_mag/sqrt(Neff),0*vel_clrs(k,:),0.3);
        for k = 1:3
            hprof(k) = plot(axp,y,vel_prof(:,k),'o-','color',vel_clrs(k,:),'markerfacecolor',vel_clrs(k,:),'linewidth',lw,'markersize',5);
        end
        hprof(4) = plot(axp,y,mag_comps,'o-','color',0.5*[1 1 1],'markerfacecolor',0.5*[1 1 1],'linewidth',lw,'markersize',5);
        hprof(5) = plot(axp,y,vel_mag,'o-','color','k','markerfacecolor','k','linewidth',lw,'markersize',5);
        xlim(axp,[0 0.525])
%         ylim(axp,prof_ylims(i,:))
        ylim(axp,[-.05 .165])
        box(axp,'on')
        xlabel(axp,'y [m]','fontsize',fs)
        ylabel(axp,'vel [m/s]')
        if i == 3
            lgd = legend(axp,hprof,{'$\textsf{u}$','$\textsf{v}$','$\textsf{w}$','$\sqrt{\textsf{u}^2+\textsf{v}^2+\textsf{w}^2}$','$|\vec{\textsf{u}}|$'},'fontsize',fs-2,'orientation','vertical','location','northeast','interpreter','latex');
            lgd.Position = lgd.Position + [0 .04 0 0];
        end

        % save figure(s)
        a = 5;
%         print(fig,sprintf('segment_summary_%02d_%02d.png',dep_nums(i),j),'-dpng','-r350')
%         print(fig_prof,sprintf('vel_prof_%02d_%02d.png',dep_nums(i),j),'-dpng','-r350')
    end
end

%%% subfunctions %%%
function hp = plotProfConf(ax,x,y_l,y_h,clr,alpha)
    % function to create shaded confidence region, inputs should be column
    % vectors
    x_patch = [x; flip(x)];
    y_patch = [y_l; flip(y_h)];
    hp = patch(ax,x_patch,y_patch,clr,'facealpha',alpha,'edgecolor','none');
end
