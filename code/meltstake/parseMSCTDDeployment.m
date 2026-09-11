function ctd = parseMSCTDDeployment(depname,varargin)
% Function to parse Aanderaa CTD data from a meltstake deployment.
%
% ctd = parseMSCTDDeployment(depname)
% ctd = parseADVDeployment(depname,t1,t2)
% ctd = parseMSCTDDeployment(depname,__,'data_folder',data_folder)
%
%
% Input
%   depname: path to the folder containing the subfolder (usually named
%       with a time stamp) which holds the raw file(s)
%   t1: (optional) start of time window to parse (datetime)
%   t2: (optional) end of time window to parse (datetime)
%   data_folder: (name-value pair) instead of searching for a data
%       directory in depname/ which holds the raw file(s), the
%       script looks directly in depname/data_folder (useful outside of meltstake
%       file tree) 
%
% Output
%   ctd: structure containing timestamps, temperature, conductivity, and
%   derived salinity
% 
% KJW
% 21 Jul 2026

%%% Check that directory exists %%%
if ~exist(depname,'dir')
    warning('directory does not exist: %s\n',depname);
    return
end

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
d = dir(depname);
if data_folder == char(1)
    data_folder = '';
    for i = 1:length(d)
        if ~(d(i).name(1)=='.') && exist(fullfile(depname,d(i).name),'dir') && exist(fullfile(depname,d(i).name,'ctd.log'),'file')
            data_folder = d(i).name;
            break
        end
    end
    if isempty(data_folder)
        error('unable to find ctd.log in %s',data_folder)
    else
        data_path = fullfile(depname,data_folder);
    end
else % use data_folder override if specified
    if exist(data_folder,'dir') && exist(fullfile(data_folder,'ctd.log'),'file')
        data_path = data_folder;
    else
        error('unable to find ctd.log in %s',data_folder)
    end
end

% find the text file in the found data folder
d2 = dir(fullfile(data_path,'raw_*.txt'));
dat_tbl = readtable(fullfile(data_path,d2(1).name));

% preallocate
ctd = struct();

% parse timestamps, C, T
ts = dat_tbl.Var1;
ctd.time = NaT(length(ts),1);
for i = 1:length(ts)
    % time
    if length(ts{i}) ~= 29 % not a full timestamp due to end of record
        break
    end
    ctd.time(i) = datetime(ts{i}(1:end-6));
end

% trim to time limits
idxt = ctd.time >= t1 & ctd.time <= t2;
ctd.time = ctd.time(idxt);

% C
ctd.C = dat_tbl.Var6(idxt);

% T
ctd.T = dat_tbl.Var8(idxt);

% calculate S
ctd.S = gsw_SP_from_C(ctd.C,ctd.T,0);

% units
C_units = dat_tbl.Var5{1}(14:end-1);
T_units = dat_tbl.Var7{1}(13:end-1);

assert(strcmp(C_units,'mS/cm'),'unexpected C units, check raw file')
assert(strcmp(T_units,'Deg.C'),'unexpected T units, check raw file')

ctd.units = {C_units, T_units, 'PSU'};