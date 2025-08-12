function rhib = parseRHIBGPS(depname,rhib_name,varargin)
% Function to parse a RHIB deployment GPS log into a convenient form.
% 
% rhib = parseRHIBGPS(depname,rhib_name)
% rhib = parseRHIBGPS(depname,rhib_name,proc_method)
% rhib = parseRHIBGPS(depname,rhib_name,t1,t2)
% rhib = parseRHIBGPS(depname,rhib_name,t1,t2,proc_method)
% 
% Input
%   depname: name of folder containing RHIB_NAME folder
%   rhib_name: name of rhib as it appears in directory
%   t1: (optional) start of time window to parse
%   t2: (optional) end of time window to parse
%   proc_method: (optional) method for calculating gps information (see
%       ROSE_code)
% 
% Output
%   rhib: structure with gps information
% 
% KJW
% 13 Sep 2022

addpath('../../../code/ROSE_code');

% Look for instrument.txt file, get dt (time coordination)
if exist(fullfile(depname,rhib_name,'instruments.txt'),'file')
    info = parseInstrumentFile(depname,rhib_name);
    fprintf('instrument file found, using time offset dt=%.1f\n',info.t_offset);
    dt = info.t_offset/86400;
    use_inst_file = true;
else
    fprintf('no instrument.txt file found\n');
    dt = 0;
    use_inst_file = false;
end

% Deal with varargin
t1 = 0;
t2 = 10^10;
proc_method = 'GPRMC groundspeed and course';
if nargin == 3 % proc_method only
    proc_method = varargin{1};
elseif nargin == 4 % time only
    t1 = varargin{1};
    t2 = varargin{2};
elseif nargin == 5 % time and proc_method
    t1 = varargin{1};
    t2 = varargin{2};
    proc_method = varargin{3};
end

% Get deployment folder
folder = extractfield(dir(fullfile(depname,rhib_name,'UBOX*')),'name');
folder = folder{1};

% Parse
fnames = extractfield(dir(fullfile(depname,rhib_name,folder,'GPS','GPS_*.log')),'name');
gps_files = fullfile(depname,rhib_name,folder,'GPS',fnames);
data = struct('gps',parse_gps(gps_files));
if ~isempty(data.gps)
    rhib = compute_vessel_vel(data,proc_method);
    rhib.time = rhib.time+dt;
    rhib.u = rhib.vel(:,1);
    rhib.v = rhib.vel(:,2);
    rhib = rmfield(rhib,'vel');
else
    rhib = struct();
    rhib.time = [];
    rhib.u = [];
    rhib.v = [];
    rhib.lon = [];
    rhib.lat = [];
end

% trim
idx = rhib.time>=t1 & rhib.time<=t2;
flds = fieldnames(rhib);
for i = 1:length(flds)
    fld = flds{i};
    rhib.(fld) = rhib.(fld)(idx);
end

% for posterity
rhib.name = rhib_name;
rhib.processing = proc_method;
rhib.description = 'Velocities are zonal/meridional in m/s. Lat/lon are decimal degrees.';
if use_inst_file
    rhib.timeCoordinationOffset = dt;
else
    rhib.timeCoordinationOffset = NaN;
end








