function hobo = parseHOBODeployment(depname,varargin)
% Function to collect HOBO files and deployment information from a
% deployment into a single structure.
% 
% hobo = parseHOBODeployment(depname)
% hobo = parseHOBODeployment(depname,t1,t2)
% hobo = parseHOBODeployment(depname,hobo_dir)
% hobo = parseHOBODeployment(depname,t1,t2,hobo_dir)
%
% Input
%   depname: name of folder containing HOBO/*.csv files
%   t1: (optional) start of time window to parse
%   t2: (optional) end of time window to parse
%   hobo_dir: (optional) directory within depname in which to look for rsk
%   files
% 
% Output
%   hobo: 1 x (number of instruments) structure
%
% KJW
% 26 Jun 2024
%
% Notes
% - added compatibility with pressure sensors, although the parsing
% algorithm could be further generalized (9/10/24)

% Parse varargin
% Optional time window
if nargin < 3
    t1 = datetime([0 1 1]);
    t2 = datetime([3000 1 1]);
else
    t1 = varargin{1};
    t2 = varargin{2};
end

% Optional sub-directory
hobo_dir = 'hobo';
if nargin == 2 || nargin == 4
    hobo_dir = varargin{end};
end

% Look for instrument.txt file
if exist(fullfile(depname,hobo_dir,'instruments.txt'),'file')
    fprintf('using instrument.txt file\n');
    hobo = parseInstrumentFile(depname,hobo_dir);
    sn = [hobo.sn];
    
else
    fprintf('no instrument.txt file found\n');
    hobo = struct();
    sn = [];
end

% Check that directory exists
if ~exist(fullfile(depname,hobo_dir),'dir')
    warning('directory does not exist: %s\n',fullfile(depname,hobo_dir));
    return
end

% Get files
fnames = {dir(fullfile(depname,hobo_dir,'*.csv')).name};
nf = length(fnames);

% Read files
for i = 1:nf
    fprintf('%d/%d\n',i,nf);
    % read header line
    fid = fopen(fullfile(depname,hobo_dir,fnames{i}),'r');
    hdr_line = fgetl(fid);
    fclose(fid);
    hdr_flds = strsplit(hdr_line,'","');
    % open file as table
    dat = readtable(fullfile(depname,hobo_dir,fnames{i}));
    % find locations of data
    for j = 1:length(hdr_flds)
        if ~isempty(regexp(hdr_flds{j},'Date Time'))
            t_col = j;
        elseif ~isempty(regexp(hdr_flds{j},'High Range')) || ~isempty(regexp(hdr_flds{j},'Low Range'))
            C_col = j;
        elseif ~isempty(regexp(hdr_flds{j},'Temp'))
            T_col = j;
        elseif ~isempty(regexp(hdr_flds{j},'Specific Conductance'))
            SC_col = j;
        elseif ~isempty(regexp(hdr_flds{j},'Salinity'))
            S_col = j;
        elseif ~isempty(regexp(hdr_flds{j},'Abs Pres'))
            P_col = j;
        end
    end
    
    % get serial number from header
    sn_start = regexp(hdr_flds{T_col},'LGR S/N: ')+9;
    sn_len = find(hdr_flds{T_col}(sn_start:end)==',',1)-1;
    sni = str2double(hdr_flds{T_col}(sn_start:sn_start+sn_len-1));

    sn_idx = find(sni==sn);    
    if isempty(sn_idx) % if serial number not found in deployment, stick it at the end and fill in metadata
        if i == 1
            sn_idx = 1;
        else
            sn_idx = length(hobo)+1;
        end
        % fill in instrument data since it's not in the instrument file
        hobo(sn_idx).sn = sni;
        hobo(sn_idx).pos = NaN;
        hobo(sn_idx).t_offset = NaN;
    end
    
    % time vector shift
    if ~isnan(hobo(sn_idx).t_offset)
        dt = hobo(sn_idx).t_offset;
    else
        dt = seconds(0);
    end
    dat{:,t_col} = dat{:,t_col}  + years(2000) + dt;

    % units
    tz = strsplit(hdr_flds{t_col},', '); % timezone
    tz = tz{2};
    hobo(sn_idx).tz = tz;
    % T
    T_unit = regexp(hdr_flds{T_col},'Temp, .* (','match');
    T_unit = T_unit{1}(7:end-2);
    
    hobo(sn_idx).units = {T_unit};

    % optional fields
    % C
    if exist('C_col','var')
        C_unit = regexp(hdr_flds{C_col},'Range, .* (','match');
        C_unit = C_unit{1}(8:end-2);
        hobo(sn_idx).units = cat(2,hobo(sn_idx).units, C_unit);
    end
    % SC
    if exist('SC_col','var')
        SC_unit = regexp(hdr_flds{SC_col},'Conductance, .* (','match');
        SC_unit = SC_unit{1}(14:end-2);
        hobo(sn_idx).units = cat(2,hobo(sn_idx).units, SC_unit);
    end
    % S
    if exist('S_col','var')
        S_unit = regexp(hdr_flds{S_col},'Salinity, .* (','match');
        S_unit = S_unit{1}(11:end-2);
        hobo(sn_idx).units = cat(2,hobo(sn_idx).units, S_unit);
    end
    % P
    if exist('P_col','var')
        P_unit = regexp(hdr_flds{P_col},'Abs Pres, .* (','match');
        P_unit = P_unit{1}(11:end-2);
        hobo(sn_idx).units = cat(2,hobo(sn_idx).units, P_unit);
    end

    % put data in output structure
    idxt = dat{:,t_col} >= t1 & dat{:,t_col} < t2;
    hobo(sn_idx).time = dat{idxt,t_col};
    hobo(sn_idx).T = dat{idxt,T_col};
    if exist('C_col','var')
        hobo(sn_idx).C = dat{idxt,C_col};
    end
    if exist('SC_col','var')
        hobo(sn_idx).SC = dat{idxt,SC_col};
    end
    if exist('S_col','var')
        hobo(sn_idx).S = dat{idxt,S_col};
    end
    if exist('P_col','var')
        hobo(sn_idx).P = dat{idxt,P_col};
    end

end
    