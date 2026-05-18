function sonar = parse881aDeployment(depname,varargin)
% Function to parse data in a Meltstake sonar881a deployment folder. The
% function relies on the files "RunData.csv" (output from binary_convert),
% "RunIndex.csv", and "configuration.json". If working outside of the
% normal Meltstake directory structure, a data folder containing the above
% files can be directly specified.
%
% sonar = parse881a(depname)
% sonar = parse881a(depname,t1,t2)
% sonar = parse881a(depname,__,'data_folder',data_folder)
%
% Input
%   depname: path to the folder containing the subfolder (usually named
%       with a time stamp) which holds the csv files
%   t1: (optional) start of time window to parse
%   t2: (optional) end of time window to parse
%   data_folder: (name-value pair) instead of searching for a data
%       directory in depname/sonar881a/ which holds the csv files, the
%       script looks directly in data_folder (useful outside of meltstake
%       file tree) 
%
% Output
%   sonar: structure containing timestamps, angle, range, and amplitudes 
%       for each of two scans completed per timestamp
%
% KJW
% 21 Dec 2025
%
% Updated and put on github 6 Feb 2026

% parse inputs
chk_data_dir = @(x) ischar(x) || isstring(x);
p = inputParser;
addOptional(p,'t1',datetime(0,1,1),@isdatetime);
addOptional(p,'t2',datetime('now'),@isdatetime);
addParameter(p,'data_folder',char(1),chk_data_dir);
parse(p,varargin{:});
t1 = p.Results.t1; t2 = p.Results.t2;
data_folder = p.Results.data_folder;

% find folder (if not specified in input)
d = dir(fullfile(depname,'sonar881a'));
if data_folder == char(1)
    data_folder = '';
    for i = 1:length(d)
        if ~(d(i).name(1)=='.') && exist(fullfile(depname,'sonar881a',d(i).name),'dir') ...
                && exist(fullfile(depname,'sonar881a',d(i).name,'RunIndex.csv'),'file') ...
                && exist(fullfile(depname,'sonar881a',d(i).name,'RunData.csv'),'file') ...
                && exist(fullfile(depname,'sonar881a',d(i).name,'configuration.json'),'file')
            data_folder = d(i).name;
            break
        end
    end
    if isempty(data_folder)
        error('unable to find RunIndex.csv, RunData.csv, and configuration.json in %s',data_folder)
    else
        data_path = fullfile(depname,'sonar881a',data_folder);
    end
else % use data_folder override if specified
    if exist(data_folder,'dir') ...
                && exist(fullfile(data_folder,'RunIndex.csv'),'file') ...
                && exist(fullfile(data_folder,'RunData.csv'),'file') ...
                && exist(fullfile(data_folder,'configuration.json'),'file')
        data_path = data_folder;
    else
        error('unable to find RunIndex.csv, RunData.csv, and configuration.json in %s',data_folder)
    end
end

dat_file = fullfile(data_path,'RunData.csv');
idx_file = fullfile(data_path,'RunIndex.csv');
cfg_file = fullfile(data_path,'configuration.json');

fid = fopen(dat_file,'r');
dat_hdrs = strsplit(fgetl(fid),',');
fclose(fid);
nh = length(dat_hdrs);

% open files
dat_tbl = readtable(dat_file);
idx_tbl = readtable(idx_file,'readvariablenames',false);
cfg = readstruct(cfg_file);

% max range
max_rng = cfg.scan.range;

% time
time = idx_tbl{:,1};
idx_time = time>=t1 & time<=t2;
time = time(idx_time);

% scan numbers
ns = size(idx_tbl,1);
nr = size(dat_tbl,1);
scan_num = nan(ns,1);
for i = 1:ns
    scan_num(i) = str2double(idx_tbl{i,3}{1}(10:end-4));
end
scan_num = scan_num(idx_time);
ns = sum(idx_time);

row_scan = nan(nr,1);
for i = 1:nr
    row_scan(i) = str2double(dat_tbl{i,1}{1}(10:end-4));
end

% number of range/angle bins
a_col = strcmp(dat_hdrs,'headposition');
dir_col = strcmp(dat_hdrs,'stepdirection');
nb = size(dat_tbl,2)-nh+1;
na = 2*length(unique(dat_tbl{:,a_col}));
scan_data = nan(ns,na,nb);
angle_data = nan(ns,na);

% read data
all_scans = dat_tbl{:,nh:end};
all_angles = dat_tbl{:,a_col};
all_dirs = dat_tbl{:,dir_col};

% all_angles(all_angles<0) = all_angles(all_angles<0)+360;
for i = 1:ns
    idx = row_scan==scan_num(i);
    scani = all_scans(idx,:);
    anglei = all_angles(idx);

    % first angle
    a0 = anglei(1);
    fidx = find(anglei==a0);
    % check that we're not in a partial last scan
    if i==ns && length(fidx)<3
        idx_final = i-1;
        scan_data = scan_data(1:idx_final,:,:);
        angle_data = angle_data(1:idx_final,:);
        time = time(1:idx_final);
        break
    end
    f1 = fidx(1);f2 = fidx(2); f3 = fidx(3);

    % flip
    if all(strcmp(all_dirs(idx),'cw'))
        flip_idx = [f2-2:-1:f1 f2-1 f3-2:-1:f2 f3-1];
        scani = scani(flip_idx,:);
        anglei = anglei(flip_idx);
    % no flip
    else
        scani = scani(f1:f3-1,:);
        anglei = anglei(f1:f3-1);
    end
    scan_data(i,:,:) = scani;
    angle_data(i,:) = anglei;
end

% unified angle axes
assert(all(all(angle_data==repmat(angle_data(1,:),[size(angle_data,1) 1]))),'angle axes do not align between scans')
a0 = angle_data(1,1);
fidx = find(angle_data(1,:)==a0);
angle_axis = angle_data(1,fidx(1):fidx(2)-1);
assert(all(angle_data(1,:)==[angle_axis angle_axis]),'angle axes do not align within scans')

% collect data into structure
sonar = struct('time',time,'angle',angle_axis','range',max_rng*linspace(0,1,nb)',...
    'scan1',scan_data(:,fidx(1):fidx(2)-1,:),'scan2',scan_data(:,fidx(2):end,:));

% do one more silly little indexing change to get the "seam" to point away
% from the ice
[~,i_jump] = max(abs(diff(angle_axis)));
idx = [i_jump:length(angle_axis) 1:i_jump-1];
sonar.angle = sonar.angle(idx);
sonar.scan1 = sonar.scan1(:,idx,:);
sonar.scan2 = sonar.scan2(:,idx,:);

