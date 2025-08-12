% Script to calculate effective DOFs/decorrelation timescale for meltstake
% segments. Apparently I didn't save any code to do this the first time.
%
% KJW
% 28 Jan 2025

clear

% load data
ms_tbl = loadMSInfo(26:28,'manualwindows');

% data location
proc_dir = 'F:meltstake/data/proc';

% collect folder names and number of segments in each deployment
folders = unique(ms_tbl.Folder,'stable');
n_segs = nan(length(folders),1);
for i = 1:length(folders)
    n_segs(i) = sum(strcmp(ms_tbl.Folder,folders{i}));
end

% set autocorrelation threshold
R0 = 1/2;

% preallocate decorrelation timescales
tau0_U = nan(sum(n_segs),1);
tau0_u = nan(sum(n_segs),1);
tau0_w = nan(sum(n_segs),1);
tau0 = nan(sum(n_segs),1);

N = nan(sum(n_segs),1);
N_star = nan(sum(n_segs),1);

u_std = nan(sum(n_segs),1);
w_std = nan(sum(n_segs),1);

% loop through deployment folders
seg = 1; % keep track of segment number
for i = 1:length(folders)
    dep_dir = fullfile(proc_dir,folders{i});

    % loop through segments
    for j = 1:n_segs(i)


        % load the stuff
        load(fullfile(dep_dir,sprintf('T%d.mat',j)))
        load(fullfile(dep_dir,sprintf('adcp%d.mat',j)))
        
        % compute timeseries to perform autocorrelations upon
        T_all = mean(T.T,2,'omitnan');
        U = vecnorm(squeeze(mean(adcp.burst.vel_ice(:,:,1:3),2,'omitnan')),2,2);
        u = squeeze(mean(adcp.burst.vel_ice(:,:,1),2,'omitnan'));
        w = squeeze(mean(adcp.burst.vel_ice(:,:,2),2,'omitnan'));
        N(seg) = length(U);

        % clean up noise just a lil
        hann_width = 0.5; % s
        U = hannFilter(U,hann_width*adcp.burst.samplerate);
        u = hannFilter(u,hann_width*adcp.burst.samplerate);
        w = hannFilter(w,hann_width*adcp.burst.samplerate);
        U_comp = u+1i*w;
        
        u_std(seg) = std(u,'omitnan');
        w_std(seg) = std(w,'omitnan');

        % autocorrelations
        max_lag = 60; % s
        kT = max_lag/seconds(mode(diff(T.time)));
        kU = max_lag*adcp.burst.samplerate;
%         kU = N(seg)-1;
        tauT = (-kT:kT)'*max_lag/max(kT);
        tauU = (-kU:kU)'/adcp.burst.samplerate;

%         R_T = laggedCrossCorr(T_all,T_all,kT);
        R_U = laggedCrossCorr(U,U,kU);
%         R_Ucomp = laggedCrossCorr(U_comp,conj(U_comp),kU);
        R_u = laggedCrossCorr(u,u,kU);
        R_w = laggedCrossCorr(w,w,kU);
        
        N_star(seg) = length(U)/sum(abs(R_U),'omitnan');
        
        tiU = findFirstZero(tauU(kU+1:end),R_U(kU+1:end)-R0);
        tiu = findFirstZero(tauU(kU+1:end),R_u(kU+1:end)-R0);
        tiw = findFirstZero(tauU(kU+1:end),R_w(kU+1:end)-R0);

        if ~isempty(tiU)
            tau0_U(seg) = tiU;
        end
        if ~isempty(tiu)
            tau0_u(seg) = tiu;
        end
        if ~isempty(tiw)
            tau0_w(seg) = tiw;
        end

        tau0(seg) = sum([u_std(seg) w_std(seg)].*[tau0_u(seg) tau0_w(seg)],'omitnan')/sum([u_std(seg) w_std(seg)].*~isnan([tau0_u(seg) tau0_w(seg)]),'omitnan');
        
        fprintf('%02d: %.2f, %.2f, %.2f (%.2f)\n',seg,tau0_U(seg),tau0_u(seg),tau0_w(seg),tau0(seg))
%         fprintf('%02d: %.1f, %5d, %.2f\n',seg,round(N_star(seg),1),N(seg),round(N_star(seg)/N(seg),2))
        
        clf(figure(1))
        hold on
        plot(tauU(kU+1:end),R_U(kU+1:end),'.-','color',colors(1))
        plot(tau0_U(seg),R0,'kx')
        plot(tauU(kU+1:end),R_u(kU+1:end),'.-','color',colors(2))
        plot(tau0_u(seg),R0,'kx')
        plot(tauU(kU+1:end),R_w(kU+1:end),'.-','color',colors(3))
        plot(tau0_w(seg),R0,'kx')
        grid on
        title(sprintf('%d',seg))

        % next segment
        seg = seg + 1;
    end
end