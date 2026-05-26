% Script to calculate candidates for integral-sacle velocity based on ADCP
% bins furthest from the wall and low-passed at various timescales.
%
% KJW
% 25 May 2026
clear

proc_dir = 'F:/meltstake/data/proc';
raw_dir = 'F:/meltstake/data/raw';

save_output = 1;

% deployment info
dep_tbl = loadMSInfo(26:28);

dep_nums = dep_tbl.Number;
n_deps = length(dep_nums);

% timescales
tau = 30.*2.^(-5:2)';

% binning
nf_u = [2 4 10 20]';
nf_u_stops = [0 80 401 1001]'/3600;

% sample volume
r_lim = [0.1 0.24];

% figure things
fs = 11;
lw = 1;
figsize = [22 24];
pad = [.05 .05 .08 .05];
shift = [0 0];
line_colors = flipud(parula(length(tau)+1));
line_colors(1,:) = [];
h_lbls = cell(length(tau)+1,1);
h_lbls{1} = 'original';
for i = 1:length(tau)
    h_lbls{i+1} = sprintf('%.1f s',round(tau(i),1));
end
pnl_lbls = {'\Psi_u [m^2 s^{-2} s^{1}]','\Psi_v [m^2 s^{-2} s^{1}]','\Psi_w [m^2 s^{-2} s^{1}]','','u [m/s]','v [m/s]','w [m/s]','(u^2+w^2)^{1/2}'};

% loop through deployments
for i = 1:n_deps
    % load
    dep_name = dep_tbl.Folder{i};
    load(fullfile(raw_dir,dep_name,'adcp','adcp.mat'));

    % update QC and transform
    if isfield(adcp.burst,'processing')
        cor_min = adcp.burst.processing.cor_min;
        amp_min = adcp.burst.processing.amp_min;
    else
        cor_min = 50;
        amp_min = 0;
    end
    adcp = msADCPTransform(adcp,cor_min,amp_min);
    
    % get rid of times not in a segment
    seg_tbl = loadMSInfo(dep_nums(i),'segments');
    n_segs = size(seg_tbl,1);
    idxt = false(size(adcp.burst.time));
    for j = 1:n_segs
        t1 = seg_tbl.Start(j);
        t2 = seg_tbl.End(j);
        idxt = idxt | (adcp.burst.time >= t1 & adcp.burst.time <= t2);
    end

    % get fully time-resolved, volume averaged velocity
    idxr = adcp.burst.range >= r_lim(1) & adcp.burst.range <= r_lim(2);
    vel = squeeze(mean(adcp.burst.vel_ice(:,idxr,[1 3 2]),2,'omitnan'));
    vel(~idxt,:,:) = nan;

    % low pass
    vel_lowpass = repmat(vel,[1 1 length(tau)]);
    for j = 1:length(tau)
        nj = 2*adcp.burst.samplerate*tau(j);
        for k = 1:size(vel,2)
            vel_lowpass(:,k,j) = hannFilter(vel_lowpass(:,k,j),nj);
        end
    end

    % psds
    vel_psd = vel(:,:,1);
    vel_psd(isnan(vel_psd)) = 0;
    [P0,f] = psd_matlab(vel_psd,adcp.burst.samplerate);
    nf = length(f);
    P = nan(nf,size(vel_lowpass,2),size(vel_lowpass,3));
    for j = 1:length(tau)
        vel_psd = vel_lowpass(:,:,j);
        vel_psd(isnan(vel_psd)) = 0;
        Pj = psd_matlab(vel_psd,adcp.burst.samplerate);
        P(:,:,j) = Pj;
    end
    % bin
    [P0bin,fbin] = progressiveBin(f,P0(:,1),nf_u,nf_u_stops);
    P0bin = repmat(P0bin,[1 3]);
    for k = 2:3
        P0bin(:,k) = progressiveBin(f,P0(:,k),nf_u,nf_u_stops);
    end
    Pbin = nan(length(fbin),size(vel_lowpass,2),size(vel_lowpass,3));
    for j = 1:length(tau)
        for k = 1:size(P,2)
            Pbin(:,k,j) = progressiveBin(f,P(:,k,j),nf_u,nf_u_stops);
        end
    end

    % set up figure
    clear fig ax h
    fig = figure(i); clf
    setFigureSize(fig,figsize);
    for j = 8:-1:1
        ax(j) = axes(fig,'position',axgridpos(4,3,j,pad,shift,'flip'));
        hold on
        box on
    end
    for j = 1:4
        ax(j+4).Position(3) = 2*ax(j+4).Position(3)+pad(1);
    end

    % plot time series
    for j = 1:3
        h(1) = plot(ax(j+4),adcp.burst.time,vel(:,j),'-','color',0.4*[1 1 1]);
        for k = 1:length(tau)
            h(k+1) = plot(ax(j+4),adcp.burst.time,vel_lowpass(:,j,k),'-','color',line_colors(k,:));
        end
    end
    plot(ax(8),adcp.burst.time,vecnorm(vel(:,[1 3]),2,2),'-','color',0.4*[1 1 1])
    for k = 1:length(tau)
        plot(ax(8),adcp.burst.time,vecnorm(vel_lowpass(:,[1 3],k),2,2),'-','color',line_colors(k,:))
    end
    linkaxes(ax(5:end),'x')

    % plot psds
    for j = 1:3
        plot(ax(j),fbin,P0bin(:,j),'-','color',0.4*[1 1 1])
        for k = 1:length(tau)
            plot(ax(j),fbin,Pbin(:,j,k),'-','color',line_colors(k,:))
        end
        set(ax(j),'xscale','log','yscale','log')
        ylim(ax(j),10.^[-8 1.5])
        xlabel(ax(j),'f [Hz]','fontsize',fs)
    end

    legend(ax(4),h,h_lbls,'fontsize',11)
    set(ax(4),'visible','off')

    pnls = addPanelLabels(ax,fs);
    for j = 1:8
        pnls(j).String = pnl_lbls{j};
        pnls(j).FontWeight = 'bold';
    end

    title(ax(5),sprintf('%s | range: [%.2f %.2f] m (%d cells)',strrep(dep_tbl.Folder{i},'_','\_'),r_lim(1),r_lim(2),sum(idxr)),'fontsize',fs+1)
    
    % save
    vel_lowpass = cat(3,vel,vel_lowpass);
    time = adcp.burst.time;
    if save_output
        save(fullfile(proc_dir,dep_tbl.Folder{i},'velocity_lowpass.mat'),'time','tau','r_lim','vel_lowpass')
    end
end
    


