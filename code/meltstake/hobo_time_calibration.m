% Script to determine hobo time offsets relative to solos
%
% 13 Oct 2025

clear

addpath('..')

% load
raw_dir = 'F:/meltstake/data/raw';
ms_tbl = loadMSInfo(28);

% hobos, solos
load(fullfile(raw_dir,ms_tbl.Folder{1},'hobo','hobo.mat'))
load(fullfile(raw_dir,ms_tbl.Folder{1},'rbr','T.mat'))

% solo sets corresponding to hobos
solo_sets = {[1 2],[2 3],[3]};

% lagged cross corr setup
k_lag = 20;
R = nan(2*k_lag+1,length(hobo));
imax = nan(1,length(hobo));

for i = 1:length(hobo)
    % trim to good data
    idx_hob = hobo(i).time>=ms_tbl.Start & hobo(i).time<=ms_tbl.End;
    idx_sol = T(solo_sets{i}(1)).time>=ms_tbl.Start & T(solo_sets{i}(1)).time<=ms_tbl.End;
    time_sol = T(solo_sets{i}(1)).time(idx_sol);
    time_hob = hobo(i).time(idx_hob);
    T_hob = hobo(i).T(idx_hob);
    
    % average solo sets
    T_sol = zeros(size(time_sol));
    for j = 1:length(solo_sets{i})
        T_sol = T_sol + T(solo_sets{i}(j)).values(idx_sol);
    end
    T_sol = T_sol/length(solo_sets{i});
    
    % upsample hobo to solo sample rate
    dt = mean(seconds(diff(time_sol)));
    T_hob_interp = interp1(time_hob,T_hob,time_sol);

    % lagged cross corr
    Ri = laggedCrossCorr(T_sol,T_hob_interp,k_lag);
    R(:,i) = Ri;
    [~,imax(i)] = max(Ri);
end
tau = dt*(-k_lag:k_lag);
T_lag = tau(imax)

% plot
figure
plot(tau,R,'.-')

