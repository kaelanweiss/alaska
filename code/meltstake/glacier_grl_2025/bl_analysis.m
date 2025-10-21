% Script to analyze boundary layer structure and time characteristics for
% each segment during the glacier meltstake deployments.
%
% KJW
% 16 Oct 2025

clear

addpath ..

proc_dir = 'F:/meltstake/data/proc';
fig_dir = 'F:/meltstake/figures/segment_summary';

load glacier_grl_2025/glacier_clrs.mat

% segment info
dep_tbl = loadMSInfo(26:28);

dep_nums = unique(dep_tbl.Number);
ndeps = length(dep_nums);

% loop through dep/segs
for i = 2:ndeps
    seg_tbl = loadMSInfo(dep_nums(i),'segments');
    nsegs = size(seg_tbl,1);
    for j = 1:nsegs
        % load adcp
        load(fullfile(proc_dir,dep_tbl.Folder{i},sprintf('adcp%d.mat',j)))
        adcp_lp = adcp;
        rmax = seg_tbl.rmax(j);
        umax = msTable2Vector(seg_tbl.u_max(j));
        y = rmax/cosd(25) - adcp.burst.range;

        % low pass beam velocities
        f_lp = 0.1*umax/(2*rmax*tand(25));
        for k = 1:adcp_lp.burst.nbeams
            for m = 1:size(adcp_lp.burst.vel,2)
                adcp_lp.burst.vel_unw(:,m,k) = lowpass(adcp_lp.burst.vel_unw(:,m,k),f_lp,adcp_lp.burst.samplerate);
            end
        end

        % perform beam transformation
        adcp_lp = msADCPTransform(adcp_lp);

        % calculate mean profiles
        vb_prof = squeeze(mean(abs(adcp.burst.vel_unw),1,'omitnan'));
        vb_prof_lp = squeeze(mean(abs(adcp_lp.burst.vel_unw),1,'omitnan'));
        vel_mag = mean(vecnorm(adcp.burst.vel_ice,2,3),'omitnan')';
        vel_mag_lp = mean(vecnorm(adcp_lp.burst.vel_ice,2,3),'omitnan')';

        fig1 = figure(1); clf
        hold on; box on

        plot(y,vb_prof,'-','linewidth',1)
        plot(y,vel_mag/2,'k-','linewidth',1)
        for k = 1:5
            plot(y,vb_prof_lp(:,k),'--','color',colors(k),'linewidth',1)
        end
        plot(y,vel_mag_lp/2,'k--','linewidth',1)
        title(sprintf('%d (%d) T_{lp}=%.1fs',dep_nums(i),j,round(1/f_lp,2)))
        a = 1;
        
        
    end
end
