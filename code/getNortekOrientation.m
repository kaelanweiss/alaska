function orientation = getNortekOrientation(dep_name,varargin)
% Function to get the pressure record from a Nortek deployment. Uses
% "Burst" or "BurstHR" records.
% 
% orientation = getNortekOrientation(dep_name)
% orientation = getNortekOrientation(dep_name,file_numbers)
% 
% Input
%   dep_name: name of deployment/folder containing /Nortek/*.mat directory
%   file_numbers: (optional) vector of file numbers to include, all files
%       included of not specified
% 
% Output
%   orientation: structure with fields
%       sn: serial number
%       t: timestamps
%       p: pressure
%       hdg: heading
%       ptc: pitch
%       rol: roll
%
% KJW
% 29 Aug 2022

% Files
d = dir(fullfile(dep_name,'Nortek','*_*.mat'));
fnames = {d.name};
nf = length(fnames);
if nf == 0
    warning('No files found in %s\n',fullfile(dep_name,'Nortek'));
    orientation = struct();
    return
end

% Get file order (dir sorts lexicographically, we want numerically)
fnums = NaN(nf,1);
for i = 1:nf
    idx1 = regexp(fnames{i},'_[0-9]*.mat')+1;
    idx2 = regexp(fnames{i},'.mat')-1;
    fnums(i) = str2double(fnames{i}(idx1:idx2));
end
[~,srt_order] = sort(fnums);
fnames = fnames(srt_order);

% select only chosen files (apply file_numbers)
if nargin == 2
    fnames = fnames(intersect([1:nf],varargin{1}));
    nf = length(fnames);
end

% Instrument config
cfg = load(fullfile(dep_name,'Nortek',fnames{1}),'Config');
cfg = cfg.Config;
if isfield(cfg,'bursthr_enable') && cfg.bursthr_enable
    prefix = 'BurstHR_';
else
    prefix = 'Burst_';
end

% Extract pressure and time
% t = cell(nf,1);
% p = cell(nf,1);
orientation = struct('sn',cfg.serialNumberDoppler,'t',[],'p',[],'hdg',[],'ptc',[],'rol',[]);
for i = 1:nf
    fprintf('%d/%d\n',i,nf);
    loadvars = load(fullfile(dep_name,'Nortek',fnames{i}),'Data');
    orientation.t = [orientation.t; loadvars.Data.([prefix 'MatlabTimeStamp'])];
    orientation.p = [orientation.p; loadvars.Data.([prefix 'Pressure'])];
    orientation.hdg = [orientation.hdg; loadvars.Data.([prefix 'Heading'])];
    orientation.ptc = [orientation.ptc; loadvars.Data.([prefix 'Pitch'])];
    orientation.rol = [orientation.rol; loadvars.Data.([prefix 'Roll'])];
%     t{i} = loadvars.Data.([prefix 'MatlabTimeStamp']);
%     p{i} = loadvars.Data.([prefix 'Pressure']);
end