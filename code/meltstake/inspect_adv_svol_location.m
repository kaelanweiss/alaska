% Script to inspect ADV distance and velocity to determine when the
% sampling volume is clear of the ice.
%
% KJW
% 18 May 2026

clear

raw_dir = 'F:/meltstake/data/raw';
proc_dir = 'F:/meltstake/data/proc';
fig_dir = 'F:/meltstake/figures/adv_sample_volume_location';

% segment info
ms_tbl = loadMSInfo(26:28,'segments');

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
for i = 1:ndeps
    nsegs = sum(ms_tbl.Number==dep_nums(i));
    dep_name = ms_tbl.Folder{ms_tbl.Number==dep_nums(i) & ms_tbl.Window==1};
    seg_start = ms_tbl.Start(ms_tbl.Number==dep_nums(i));
    seg_end = ms_tbl.End(ms_tbl.Number==dep_nums(i));
    
    % load data
    load(fullfile(raw_dir,dep_name,'adv','adv.mat'))
    load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'))
    load(fullfile(proc_dir,dep_name,'adv','pck_edges.mat'))

    % perform ADV transform
    adv = msADVTransform(adv,adcp.attitude);
    adv = burstAverageADV(adv);

    % normalize backscatter by mean profile
    % adv.pck_amp = adv.pck_amp./repmat(mean_prof',[size(adv.pck_amp,1) 1 3]);

    % get some backscatter values for clims
    pck_amp_min = mean(adv.pck_amp(:,adv.pck_dist(2,:) > 30 & adv.pck_dist(2,:) < 50,:),'all','omitnan');

    % instrument velocity components
    fig2 = figure(2); clf
    setFigureSize(fig2,figsize.*[1 3/4]);
    clear ax2
    for j = 1:3
        ax2(j) = axes(fig2,'position',axgridpos(3,1,j,pad,shift));
        hold on
        box on
        % plot(ax2(j),adv.time,adv.vel_ice,'-','color',0.5*[1 1 1])
        plot(ax2(j),extrema(adv.time),[0 0],'k-')
        plot(ax2(j),adv.time,adv.vel_ice(:,j),'-','color',vel_clrs(j,:),'linewidth',lw)
    end
    linkaxes(ax2,'x')

    % backscatter amplitudes
    fig4 = figure(4); clf
    setFigureSize(fig4,[27 9.5]);
    clear ax4
    for j = 1:3
        ax4(j) = axes(fig4,'position',axgridpos(1,3,j,[.03 .05 .1 .1],[0 0]));
        hold on
        pcolor(adv.pck_time,adv.pck_dist(2,:),adv.pck_amp(:,:,j)')
        shading flat
        % cmocean('rain')
        clim([pck_amp_min 0.9*max(adv.pck_amp,[],'all')])

        % ice edges
        idx_in = edges.pos(:,j) < (157 + adv.sampling_volume/2);
        plot(edges.time(~idx_in),edges.pos(~idx_in,j),'.','color',0*[1 1 1])
        plot(edges.time(idx_in),edges.pos(idx_in,j),'r.')

        % sample volume
        plot(extrema(adv.time),157*[1 1],'k-','linewidth',lw)
        plot(extrema(adv.time),157*[1 1]-adv.sampling_volume/2,'k--','linewidth',lw-0.2)
        plot(extrema(adv.time),157*[1 1]+adv.sampling_volume/2,'k--','linewidth',lw-0.2)
    end
    ylabel(ax4(1),'distance [mm]','fontsize',fs)
    linkaxes(ax4)
    ylim(ax4(1),[50 300])
    
    lbls = addPanelLabels(ax4,fs);
    for j = 1:3
        lbls(j).String = sprintf('beam %d',j);
    end
    linkaxes([ax2 ax4],'x')

    for j = 1:nsegs
        xlim(ax2(1),[seg_start(j) seg_end(j)])
        title(ax4(2),sprintf('%s | %d(%d)',strrep(dep_name,'_','\_'),dep_nums(i),j))

        % save figures
        print(fig4,fullfile(fig_dir,sprintf('adv_smpvol_location_%02d_%02d.png',dep_nums(i),j)),'-dpng','-r300')
    end
    a=1;
end
