% Script to plot processed data from each meltstake segment for easy
% inspection -- velocity profiles.
%
% KJW
% 9 Oct 2025
clear

proc_dir = 'F:/meltstake/data/proc';
fig_dir = 'F:/meltstake/figures/segment_summary';

load glacier_grl_2025/glacier_clrs.mat

% segment info
ms_tbl = loadMSInfo(26:28,'segments');

dep_nums = unique(ms_tbl.Number);
ndeps = length(dep_nums);

% plotting
lw = 1.5;
ms = 4;
fs = 12;
figsize = [13 8];
pad = [.05 .02 .15 .07];
shift = [0 0];
vel_comp_order = [1 3 2];
vel_lbls = {'u','v','w','U_{max}'};

% loop through deployment
for i = 1:ndeps
    nsegs = sum(ms_tbl.Number==dep_nums(i));
    dep_name = ms_tbl.Folder{ms_tbl.Number==dep_nums(i) & ms_tbl.Window==1};
    for j = 1:nsegs
        % load data
        load(fullfile(proc_dir,dep_name,sprintf('adcp%d.mat',j)))

        % compute profile things
        rmax = ms_tbl.rmax(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        y = rmax/cosd(25) - adcp.burst.range;
        vel_prof = squeeze(mean(adcp.burst.vel_ice(:,:,vel_comp_order),1,'omitnan')); % component means
        mag_comps = vecnorm(vel_prof,2,2); % magnitude of mean components
        tau = ms_tbl.tau_decorr(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j); % decorrelation
        Neff = 60*ms_tbl.Duration(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j)/ms_tbl.tau_decorr(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        if isnan(Neff) % failsafe, it's wrong but whatever
            Neff = length(adcp.burst.time);
        end
        prof_std = squeeze(std(adcp.burst.vel_ice(:,:,vel_comp_order),0,1,'omitnan'))/sqrt(Neff); % component standard error
        dmag_comps = sqrt(sum((vel_prof.*prof_std).^2,2))./mag_comps;
        vel_mag = mean(vecnorm(adcp.burst.vel_ice,2,3),'omitnan')';
        std_mag = std(vecnorm(adcp.burst.vel_ice,2,3),'omitnan')';

        % plot
        fig = figure(1); clf
        setFigureSize(fig,figsize);
        ax = axes(fig,'position',axgridpos(1,1,1,[0 0 .15 .15],[-.03 0]),'fontsize',fs-2);
        hold on; box on
        clear h
        % plot zero line
        plot(ax,[0 1],[0 0],'--','color',0.*[1 1 1],'linewidth',lw)
        % plot confidence regions
        for k = 1:3
            plotConfidenceRegion(ax,y,vel_prof(:,k)-prof_std(:,k),vel_prof(:,k)+prof_std(:,k),vel_clrs(k,:),.3);
        end
        plotConfidenceRegion(ax,y,mag_comps-dmag_comps,mag_comps+dmag_comps,0.5*[1 1 1],.3);
        plotConfidenceRegion(ax,y,vel_mag-std_mag/sqrt(Neff),vel_mag+std_mag/sqrt(Neff),0*[1 1 1],.3);
        % plot profiles
        for k = 1:3
            h(k) = plot(ax,y,vel_prof(:,k),'o-','color',vel_clrs(k,:),'markerfacecolor',vel_clrs(k,:),'linewidth',lw,'markersize',ms);
        end
        h(4) = plot(ax,y,mag_comps,'o-','color',0.5*[1 1 1],'markerfacecolor',0.5*[1 1 1],'linewidth',lw,'markersize',ms);
        h(5) = plot(ax,y,vel_mag,'o-','color','k','markerfacecolor','k','linewidth',lw,'markersize',ms);
        % labels and things
        v_min = min(vel_prof-prof_std,[],'all');
        v_max = max(vel_mag+std_mag/sqrt(Neff));
%         ylim(ax,[floor(v_min*20)/20 ceil(v_max*20)/20])
        xlim(ax,[0 1.05*max(y)])
        xlabel(ax,'y [m]','fontsize',fs)
        ylabel(ax,'vel [m/s]')
        lgd = legend(ax,h,{'$\textsf{u}$','$\textsf{v}$','$\textsf{w}$','$\sqrt{\textsf{u}_i^2}$','$|\vec{\textsf{u}}|$'},'fontsize',fs-2,'orientation','vertical','interpreter','latex');
        lgd.Position(1:2) = [sum(ax.Position([1 3]))+.01 ax.Position(2)+.1];
        % title
        depth = ms_tbl.depth(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        dur = ms_tbl.Duration(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j);
        m_obs = msTable2Vector(ms_tbl.m(ms_tbl.Number==dep_nums(i) & ms_tbl.Window==j))*.24;
        title(ax(1),sprintf('dep %d seg %d | %.1f m, %d min, %.1f m/day',dep_nums(i),j,round(depth,1),round(dur),round(m_obs,1)),'fontsize',fs)
        
        % see if attenuation of velocity magnitude is due to noncoherent
        % motion near the wall
        fig2 = figure(2); clf
        ax2 = axes(fig2);
        hold on; box on
        for k = 1:5
            plot(ax2,rmax/cosd(25)-adcp.burst.range,squeeze(mean(abs(adcp.burst.vel_unw(:,:,k)),1,'omitnan')),'-','color',colors(k))
        end
        plot(ax2,y,vel_mag,'o-','color','k','markerfacecolor','k','linewidth',lw,'markersize',ms);
        
        % save figure
        a = 1;
%         print(fig,fullfile(fig_dir,sprintf('velocity_profiles_%02d_%02d.png',dep_nums(i),j)),'-dpng','-r300')
    end
end
