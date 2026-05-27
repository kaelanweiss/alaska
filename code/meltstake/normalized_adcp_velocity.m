% Script to make pcolor plots of ADCP velocity normalized by integral scale
% velocity created by "integral_scale_velocity.m" script.
%
% KJW
% 26 May 2026
clear

proc_dir = 'F:/meltstake/data/proc';
raw_dir = 'F:/meltstake/data/raw';
fig_dir = 'F:/meltstake/figures';

% deployment info
dep_tbl = loadMSInfo(26:28);

dep_nums = dep_tbl.Number;
n_deps = length(dep_nums);

% loop through deployments
for i = 1:n_deps
    int_vel = load(fullfile(proc_dir,dep_tbl.Folder{i},'velocity_lowpass.mat'));
    int_vel.tau = [0; int_vel.tau];

    % loop through 