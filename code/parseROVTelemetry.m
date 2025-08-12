function rov = parseROVTelemetry(depname,varargin)
% Function to parse ROV telemetry .csv files.
% 
% rov = parseROVTelemetry(depname)
% rov = parseROVTelemetry(depname,flds_in)
% rov = parseROVTelemetry(depname,flds_in,flds_out)
% rov = parseROVTelemetry(depname,...,t1,t2)
%
% Input
%   depname: path to folder containing "QGroundControl/*.csv"
%   flds_in: (optional) additional fields to parse from telemetry log, if
%       "flds_out" not specified, field names are copied to the output
%       structure
%   flds_out: (optional) field names corresponding to "flds_in" in the
%       output structure
%   t1: (optional) start of time window to parse as datenum
%   t2: (optional) end of time window to parse as datenum
% 
% Output
%   rov: structure of parsed telemetry
% 
% KJW
% 29 Aug 2022

% Fields to parse and their output names
FLDS_IN = {'Timestamp',...
           'roll',...
           'pitch',...
           'heading',...
           'climbRate',...
           'altitudeRelative',...
           'altitudeAMSL',...
           'throttlePct'};
       
FLDS_OUT = FLDS_IN;
FLDS_OUT{1} = 'time';

t1 = 0;
t2 = datenum(3000,1,1);

% Optional: add other fields
if nargin == 2 % only input fields, no time frame
    FLDS_IN = cat(2,FLDS_IN,varargin{1});
    FLDS_OUT = cat(2,FLDS_OUT,varargin{1});
elseif nargin == 3 && iscell(varargin{1}) && iscell(varargin{2}) % input and output fields, no time frame
    FLDS_IN = cat(2,FLDS_IN,varargin{1});
    FLDS_OUT = cat(2,FLDS_OUT,varargin{2});
elseif nargin == 3 && isnumeric(varargin{1}) && isnumeric(varargin{2})  % no input/output fields, time frame
    t1 = varargin{1};
    t2 = varargin{2};
elseif nargin == 4 % input and time frame
    FLDS_IN = cat(2,FLDS_IN,varargin{1});
    FLDS_OUT = cat(2,FLDS_OUT,varargin{1});
    t1 = varargin{2};
    t2 = varargin{3};
elseif nargin == 5 % input, output, and time frame
    FLDS_IN = cat(2,FLDS_IN,varargin{1});
    FLDS_OUT = cat(2,FLDS_OUT,varargin{2});
    t1 = varargin{3};
    t2 = varargin{4};
end
nflds = length(FLDS_IN);

rov = struct();

% Check that directory exists
if ~exist(fullfile(depname,'QGroundControl'),'dir')
    warning('directory does not exist: %s\n',fullfile(depname,'QGroundControl'));
    return
end

% Look for instrument.txt file
if exist(fullfile(depname,'QGroundControl','instruments.txt'),'file')
    info = parseInstrumentFile(depname,'QGroundControl');
    fprintf('instrument file found, using time offset dt=%.1f\n',info.t_offset);
    dt = info.t_offset/86400;
    rov.timeCoordinationOffset = dt;
else
    fprintf('instrument file not found\n');
    dt = 0;
    rov.timeCoordinationOffset = NaN;
end

% Get files to parse
d = dir(fullfile(depname,'QGroundControl','*.csv'));
fnames = {d.name};
nf = length(fnames);

% Count number of lines
nlines = -nf;
for i = 1:nf
    fid = fopen(fullfile(depname,'QGroundControl',fnames{i}),'r');
    while ~feof(fid)
        line = fgetl(fid);
        try
            dn = readTimeStamp(line);
        catch
            %fprintf('%s\n',line);
            dn = 0;
        end
        if dn>=t1 && dn<=t2
            nlines = nlines + 1;
        end
    end
    fclose(fid);
end
    

% Instantiate output
for i = 1:nflds
    rov.(FLDS_OUT{i}) = NaN(nlines,1);
end

% Parse
pos = 1;
for i = 1:nf
    fid = fopen(fullfile(depname,'QGroundControl',fnames{i}),'r');
    fld_lbls = strsplit(fgetl(fid),',');
    
    % attempt to find positions of wanted fields
    fld_idx = NaN(1,nflds);
    for j = 1:nflds
        srch = find(strcmp(FLDS_IN{j},fld_lbls));
        if length(srch)==1
            fld_idx(j) = srch;
        end
    end
    
    % loop through lines
    while ~feof(fid)
        if ~mod(pos,100)
            fprintf('%d/%d\n',pos,nlines);
        end
        line = fgetl(fid);
        dn = readTimeStamp(line)+dt;
        if dn>=t1 && dn<=t2
            vals = strsplit(line,',');
            rov.(FLDS_OUT{1})(pos) = dn;
            for j = 2:nflds
                if ~isnan(fld_idx(j))
                    valj = str2double(vals(fld_idx(j)));
                    rov.(FLDS_OUT{j})(pos) = valj;
                end
            end
            pos = pos + 1;
        end
    end
    
    fclose(fid);
end
fprintf('%d/%d\n',pos-1,nlines);

% calculated fields
rov.yawRate = [medianFilter(diff(rov.heading)./round(86400*diff(rov.time),1),3); NaN];
rov.pitchRate = [medianFilter(diff(rov.pitch)./round(86400*diff(rov.time),1),3); NaN];
rov.rollRate = [medianFilter(diff(rov.roll)./round(86400*diff(rov.time),1),3); NaN];

end

function dn = readTimeStamp(line)
    ts = line(1:find(line==',',1)-1);
    dn = datenum(ts,'yyyy-mm-dd HH:MM:SS.fff');
end

