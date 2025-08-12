% Script to plot comprehensive ADV data from each meltstake segment for easy
% inspection.
%
% KJW
% 29 Jul 2025

clear

raw_dir = 'F:/meltstake/data/raw';
fig_dir = 'F:/meltstake/figures/adv_velocity';

% segment info
ms_tbl = loadMSInfo(26:28,'manualwindows');

dep_nums = unique(ms_tbl.Number);
ndeps = length(dep_nums);

% plotting
lw = 1;
fs = 12;
figsize = [15 18];
pad = [.05 .02 .1 .07];
shift = [0 0];
qc_flds = {'vel','cor','snr','amp'};
vel_comp_order = [1 3 2];
vel_lbls = {'u','v (x2)','w','U_{max}'};
load egu_2025\clrs_egu.mat
load F:/adv/mean_profile.mat

% loop through deployment
for i = 3:ndeps
    nsegs = sum(ms_tbl.Number==dep_nums(i));
    dep_name = ms_tbl.Folder{ms_tbl.Number==dep_nums(i) & ms_tbl.Window==1};
    seg_start = ms_tbl.Start(ms_tbl.Number==dep_nums(i));
    seg_end = ms_tbl.End(ms_tbl.Number==dep_nums(i));
    
    % load data
    load(fullfile(raw_dir,dep_name,'adv','adv.mat'))
    load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'))

    % perform ADV transform
    adv = msADVTransform(adv,adcp.attitude);
    adv = burstAverageADV(adv);

    % normalize backscatter by mean profile
    adv.pck_amp = adv.pck_amp./repmat(mean_prof',[size(adv.pck_amp,1) 1 3]);

    % get some backscatter values
    pck_amp_min = mean(adv.pck_amp(:,adv.pck_dist(2,:) > 30 & adv.pck_dist(2,:) < 50,:),'all','omitnan');

    % QC figure
    fig1 = figure(1); clf
    setFigureSize(fig1,figsize);
    clear ax1
    for j = 1:length(qc_flds)
        ax1(j) = axes(fig1,'position',axgridpos(length(qc_flds),1,j,pad,shift));
        plot(ax1(j),adv.time,adv.(qc_flds{j}),'-')
    end
    linkaxes(ax1,'x')

    % instrument velocity components
    fig2 = figure(2); clf
    setFigureSize(fig2,figsize.*[1 3/4]);
    clear ax2
    for j = 1:3
        ax2(j) = axes(fig2,'position',axgridpos(3,1,j,pad,shift));
        hold on
        box on
        plot(ax2(j),adv.time,adv.vel_ice,'-','color',0.5*[1 1 1])
        plot(ax2(j),extrema(adv.time),[0 0],'k-')
        plot(ax2(j),adv.time,adv.vel_ice(:,j),'-','color',vel_clrs(j,:),'linewidth',lw)
    end
    linkaxes(ax2,'x')

    % burst-average velocity
    fig3 = figure(3); clf
    setFigureSize(fig3,figsize.*[1 3/4]);
    clear ax3
    for j = 1:3
        ax3(j) = axes(fig3,'position',axgridpos(3,1,j,pad,shift));
        hold on
        box on
        for k = setxor(1:3,j)
            plotConfidenceRegion(ax3(j),adv.time_avg,adv.vel_avg(:,k)-adv.vel_std(:,k),adv.vel_avg(:,k)+adv.vel_std(:,k),0.5*[1 1 1],0.25)
        end
        plotConfidenceRegion(ax3(j),adv.time_avg,adv.vel_avg(:,j)-adv.vel_std(:,j),adv.vel_avg(:,j)+adv.vel_std(:,j),vel_clrs(j,:),0.65)
    end
    for j = 1:3
        plot(ax3(j),adv.time_avg,adv.vel_avg,'-','color',0.5*[1 1 1])
        plot(ax3(j),extrema(adv.time),[0 0],'k-')
        plot(ax3(j),adv.time_avg,adv.vel_avg(:,j),'-','color',vel_clrs(j,:),'linewidth',lw)
    end
    linkaxes(ax3,'x')

    % backscatter amplitudes
    fig4 = figure(4); clf
    setFigureSize(fig4,[27 9.5]);
    clear ax4
    for j = 1:3
        ax4(j) = axes(fig4,'position',axgridpos(1,3,j,[.02 .05 .1 .1],[0 0]));
        hold on
        pcolor(adv.pck_time,adv.pck_dist(2,:),adv.pck_amp(:,:,j)')
        shading flat
        cmocean('mat')
        clim([pck_amp_min 0.9*max(adv.pck_amp,[],'all')])
        plot(extrema(adv.time),157*[1 1],'r-')
        plot(extrema(adv.time),157*[1 1]-adv.sampling_volume/2,'r--')
        plot(extrema(adv.time),157*[1 1]+adv.sampling_volume/2,'r--')
    end
    linkaxes(ax4)
    ylim(ax4(1),[50 250])
    
    linkaxes([ax1 ax2 ax3 ax4],'x')

    for j = 1:nsegs
        xlim(ax1(1),[seg_start(j) seg_end(j)])
        a = 1;
        % save figures
        %print(fig,fullfile(fig_dir,sprintf('segment_summary_%02d_%02d.png',dep_nums(i),j)),'-dpng','-r300')
    end
end
