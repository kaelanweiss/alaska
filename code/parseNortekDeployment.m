function adcp = parseNortekDeployment(depname,varargin)
% Function to parse a Nortek Signature deployment of MAT files. A time 
% window can be specified. This should work for all Signature ADCPs in all
% configurations. No processing is performed (e.g. QC, velocity 
% transformations). The script expects the .mat files output from MIDAS to
% be located in a folder "Nortek"
%
% adcp = parseNortekDeployment(depname)
% adcp = parseNortekDeployment(depname,flds)
% adcp = parseNortekDeployment(depname,t1,t2)
% adcp = parseNortekDeployment(depname,flds,t1,t2)
% 
% Input
%   depname: path to folder containing a subfolder 'Nortek/' which contains
%       the ADCP .mat files from MIDAS, i.e. "(depname)/Nortek/*.mat"
%   flds: (optional) cell array of field names to parse, options are 'burst', 'bt', 
%       'echo', and 'none', if empty all available fields are parsed
%   t1: (optional) start of time window to parse
%   t2: (optional) end of time window to parse
% 
% Output
%   adcp: structure containing configuration information and data, fields
%       correspond to the configured sampling modes (e.g. burst, bottom
%       tracking, echosounder) as well as instrument configuration 
%       information
%
% KJW
% 1 Sep 2022
%
% NOTE: changed dir() call to include all *.mat files instead of *_*.mat
% NOTE: 30 Jan 2024 switched from datenum to datetime
%
% 12 Jan 2023

% Instantiate output
adcp = struct();

%%% Check that directory exists %%%
if exist(fullfile(depname,'Nortek'),'dir')
    subdir = 'Nortek';
elseif exist(fullfile(depname,'mat'),'dir')
    subdir = 'mat';
else
    warning('Nortek or mat directory does not exist in %s\n',depname);
    return
end

%%% Look for instrument.txt file %%%
if exist(fullfile(depname,subdir,'instruments.txt'),'file')
    info = parseInstrumentFile(depname,subdir);
    fprintf('instrument file found, using time offset dt=%.1f\n',seconds(info.t_offset));
    if ~isnan(info.t_offset)
        dt = info.t_offset;
    else
        dt = seconds(0);
    end
    use_inst_file = true;
else
    fprintf('instrument file not found\n');
    dt = seconds(0);
    use_inst_file = false;
end

%%% Files %%%
d = dir(fullfile(depname,subdir,'*.mat'));
fnames = {d.name}';
nf = length(fnames);
if nf == 0
    warning('No files found in %s\n',fullfile(depname));
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

%%% Instrument config %%%
cfg = load(fullfile(depname,subdir,fnames{1}),'Config'); cfg = cfg.Config;
adcp.cfg = struct();
adcp.cfg.sn = cfg.serialNumberDoppler;
adcp.cfg.model = cfg.instrumentName;
adcp.cfg.pressureoffset = cfg.pressureOffset;
if use_inst_file
    adcp.cfg.timeCoordinationOffset = info.t_offset;
else
    adcp.cfg.timeCoordinationOffset = NaN;
end
t1 = datetime(cfg.startCollectionTime_seconds,'convertfrom','epochtime');
t2 = datetime(ceil(cfg.endCollectionTime_seconds+cfg.endCollectionTime_subseconds),'convertfrom','epochtime');

% check if burst enabled
if isfield(cfg,'burst_enable') && cfg.burst_enable
    adcp.cfg.burst_enable = cfg.burst_enable;
else
    adcp.cfg.burst_enable = 0;
end

% check for burst HR mode
if isfield(cfg,'bursthr_enable') && cfg.bursthr_enable
    burst_prefix = 'BurstHR_';
    adcp.cfg.bursthr_enable = cfg.bursthr_enable;
else
    burst_prefix = 'Burst_';
    adcp.cfg.bursthr_enable = 0;
end

% check if bottom tracking enabled
if isfield(cfg,'bt_nBeams')
    adcp.cfg.bt_enable = 1;
    bt_prefix = 'BurstBT_';
else
    adcp.cfg.bt_enable = 0;
end

% check if echosounder enabled
if isfield(cfg,'echo_enable') && cfg.echo_enable
    adcp.cfg.echo_enable = cfg.echo_enable;
    echo_prefix = ['Echo1' num2str(cfg.echo_frequency1) '_0kHz_'];
else
    adcp.cfg.echo_enable = 0;
end

%%% time frame, parse varargin %%%
flds = {};
switch nargin
    case 1
        t1_parse = t1;
        t2_parse = t2;
    case 2   
        flds = varargin{1};
    case 3
        t1_parse = varargin{1}-dt;
        t2_parse = varargin{2}-dt;
    case 4
        flds = varargin{1};
        t1_parse = varargin{2}-dt;
        t2_parse = varargin{3}-dt;
end
if t1_parse > t1
    t1 = t1_parse;
end
if t2_parse < t2
    t2 = t2_parse;
end
duration = seconds(t2-t1);

%%% Decide what to parse %%%
maxRate = 32;
if adcp.cfg.burst_enable
    maxRate = maxRate/2;
    if cfg.burst_nBeams>4
        maxRate = maxRate/2;
    end
end
if adcp.cfg.bt_enable
    maxRate = maxRate/2;
end
if adcp.cfg.echo_enable
    maxRate = maxRate/2;
end
nt_max = duration*maxRate+100;

% fields
parse_burst = 0;
parse_bt = 0;
parse_echo = 0;
if ~isempty(flds) && all(strcmp('none',flds))
    % do nothing, keep all parse flags false
else
    % burst
    if adcp.cfg.burst_enable && (any(strcmp('burst',flds)) || isempty(flds))
        parse_burst = 1;
        adcp.burst = struct();
        adcp.burst.POS = 0;
        % config info
        adcp.burst.samplerate = cfg.burst_sampleRate;
        adcp.burst.ncells = cfg.([lower(burst_prefix) 'nCells']);
        adcp.burst.nbeams = cfg.burst_nBeams;
        adcp.burst.cellsize = cfg.([lower(burst_prefix) 'cellSize']);
        adcp.burst.blankingdist = cfg.([lower(burst_prefix) 'blankingDistance']);
        adcp.burst.coordinates = cfg.burst_coordSystem;
        adcp.burst.beam2xyz = reshape(cfg.burst_beam2xyz,[4 4])';
        adcp.burst.velrange = cfg.burst_velocityRange;
        if adcp.cfg.bursthr_enable
            adcp.burst.hrlag = cfg.bursthr_lag;
        end
        % preallocate
        % nt x ncells x nbeams
        adcp.burst.time = NaT(nt_max,1);
        adcp.burst.vel = NaN(nt_max,adcp.burst.ncells,adcp.burst.nbeams);
        adcp.burst.cor = NaN(nt_max,adcp.burst.ncells,adcp.burst.nbeams);
        adcp.burst.amp = NaN(nt_max,adcp.burst.ncells,adcp.burst.nbeams);
    end
    % bt
    if adcp.cfg.bt_enable && (any(strcmp('bt',flds)) || isempty(flds))
        parse_bt = 1;
        adcp.bt = struct();
        adcp.bt.POS = 0;
        % config info
        adcp.bt.ncells = cfg.bt_nCells;
        adcp.bt.nbeams = cfg.bt_nBeams;
        adcp.bt.cellsize = cfg.bt_cellSize;
        adcp.bt.maxrange = cfg.bt_range;
        adcp.bt.velocityrange = cfg.bt_velocityRange;
        % preallocate
        % nt x nbeams
        adcp.bt.time = NaT(nt_max,1);
        adcp.bt.vel = NaN(nt_max,adcp.bt.nbeams);
        adcp.bt.distance = NaN(nt_max,adcp.bt.nbeams);
        adcp.bt.fom = NaN(nt_max,adcp.bt.nbeams);
    end
    % echo
    if adcp.cfg.echo_enable && (any(strcmp('echo',flds)) || isempty(flds))
        parse_echo = 1;
        adcp.echo = struct();
        adcp.echo.POS = 0;
        % config info
        adcp.echo.ncells = cfg.echo_numBins;
        adcp.echo.cellsize = cfg.echo_binSize;
        adcp.echo.blankingdist = cfg.echo_blanking;
        adcp.echo.frequency = cfg.echo_frequency1;
        % preallocate
        % nt x ncells
        adcp.echo.time = NaT(nt_max,1);
        adcp.echo.amp = NaN(nt_max,adcp.echo.ncells);
    end    
end

% preallocate attitude
adcp.attitude = struct();
adcp.attitude.POS = 0;
adcp.attitude.time = NaT(nt_max,1);
adcp.attitude.pressure = NaN(nt_max,1);
adcp.attitude.heading = NaN(nt_max,1);
adcp.attitude.pitch = NaN(nt_max,1);
adcp.attitude.roll = NaN(nt_max,1);
adcp.attitude.batteryvoltage = NaN(nt_max,1);
adcp.attitude.ahrsgyro = NaN(nt_max,3);
adcp.attitude.ahrsrotationmatrix = NaN(nt_max,9);
if adcp.cfg.burst_enable
    att_prefix = burst_prefix;
elseif adcp.cfg.bt_enable
    att_prefix = bt_prefix;
elseif adcp.cfg.echo_enable
    att_prefix = echo_prefix;
end

%%% Parse %%%
t1 = datenum(t1);
t2 = datenum(t2);
for i = 1:nf
    fprintf('%d/%d\n',i,nf);
    % load MIDAS MAT file
    data = load(fullfile(depname,subdir,fnames{i}),'Data'); data = data.Data;
    if i == 1
        adcp.units = data.Units;
    end
    % 1 attitude
    in_slice = data.([att_prefix 'MatlabTimeStamp'])>=t1 & data.([att_prefix 'MatlabTimeStamp'])<=t2;
    n = length(find(in_slice));
    out_slice = [1:n]+adcp.attitude.POS;
    adcp.attitude.time(out_slice) = datetime(data.([att_prefix 'MatlabTimeStamp'])(in_slice),'convertfrom','datenum') + dt;
    adcp.attitude.pressure(out_slice) = data.([att_prefix 'Pressure'])(in_slice);
    adcp.attitude.heading(out_slice) = data.([att_prefix 'Heading'])(in_slice);
    adcp.attitude.pitch(out_slice) = data.([att_prefix 'Pitch'])(in_slice);
    adcp.attitude.roll(out_slice) = data.([att_prefix 'Roll'])(in_slice);
    adcp.attitude.batteryvoltage(out_slice) = data.([att_prefix 'Battery'])(in_slice);
    adcp.attitude.ahrsgyro(out_slice,:) = [data.([att_prefix 'AHRSGyroX'])(in_slice),...
                                           data.([att_prefix 'AHRSGyroY'])(in_slice),...
                                           data.([att_prefix 'AHRSGyroZ'])(in_slice)];
    adcp.attitude.ahrsrotationmatrix(out_slice,:) = data.([att_prefix 'AHRSRotationMatrix'])(in_slice,:);
    adcp.attitude.POS = adcp.attitude.POS + n;
    
    % 2 burst
    if parse_burst
        if i == 1
            adcp.burst.range = cast(data.([burst_prefix 'Range']),'double');
        end
        in_slice = data.([burst_prefix 'MatlabTimeStamp'])>=t1 & data.([burst_prefix 'MatlabTimeStamp'])<=t2;
        n = length(find(in_slice));
        out_slice = [1:n]+adcp.burst.POS;
        adcp.burst.time(out_slice) = datetime(data.([burst_prefix 'MatlabTimeStamp'])(in_slice),'convertfrom','datenum') + dt;
        if adcp.burst.nbeams==5 % trim 5th beam if longer
            in_slice5 = data.(['I' burst_prefix 'MatlabTimeStamp'])>=t1 & data.(['I' burst_prefix 'MatlabTimeStamp'])<=t2;
            n5 = length(find(in_slice5));
        end
        for fld_loop = {'Vel','Cor','Amp'}
            fld = fld_loop{1};
            for j = 1:adcp.burst.nbeams
                if j < 5
                    adcp.burst.(lower(fld))(out_slice,:,j) = data.([burst_prefix fld 'Beam' num2str(j)])(in_slice,:);
                else
                    tmp5 = data.(['I' burst_prefix fld 'Beam' num2str(j)])(in_slice5,:);
                    if n5 > n % make sure 5th beam data vector is same length (not guaranteed)
                        tmp5 = tmp5(1:n,:);
                    elseif n5 < n
                        tmp5 = [tmp5; nan(n-n5,size(tmp5,2))];
                    end
                    adcp.burst.(lower(fld))(out_slice,:,j) = tmp5;
                end
            end
        end
        adcp.burst.POS = adcp.burst.POS + n;
    end
    
    % 3 bt
    if parse_bt
        in_slice = data.([bt_prefix 'MatlabTimeStamp'])>=t1 & data.([bt_prefix 'MatlabTimeStamp'])<=t2;
        n = length(find(in_slice));
        out_slice = [1:n]+adcp.bt.POS;
        adcp.bt.time(out_slice) = datetime(data.([bt_prefix 'MatlabTimeStamp'])(in_slice),'convertfrom','datenum') + dt;
        for fld_loop = {'Vel','Distance','FOM'}
            fld = fld_loop{1};
            for j = 1:adcp.bt.nbeams
                adcp.bt.(lower(fld))(out_slice,j) = data.([bt_prefix fld 'Beam' num2str(j)])(in_slice);
            end
        end
        adcp.bt.POS = adcp.bt.POS + n;
    end

    % 4 echo
    if parse_echo
        if i == 1
            adcp.echo.range = data.([echo_prefix 'Range']);
        end
        in_slice = data.([echo_prefix 'MatlabTimeStamp'])>=t1 & data.([echo_prefix 'MatlabTimeStamp'])<=t2;
        n = length(find(in_slice));
        out_slice = [1:n]+adcp.echo.POS;
        adcp.echo.time(out_slice) = datetime(data.([echo_prefix 'MatlabTimeStamp'])(in_slice),'convertfrom','datenum') + dt;
        adcp.echo.amp(out_slice,:) = data.([echo_prefix 'Amp' num2str(adcp.echo.frequency) '_0kHz'])(in_slice,:);
        adcp.echo.POS = adcp.echo.POS + n;
    end
end

% trim excess preallocated space
fld_names = fieldnames(adcp);
idx_cfg_units = strcmp(fld_names,'cfg') | strcmp(fld_names,'units');
fld_names = fld_names(~idx_cfg_units)';

for fld_loop = fld_names % loop through adcp fields
    fld = fld_loop{1};
    n_trim = adcp.(fld).POS;
    flds2 = fieldnames(adcp.(fld));
    for i = 1:length(flds2) % loop through adcp subfields within each field (these are the time vectors)
        if size(adcp.(fld).(flds2{i}),1) == nt_max
            % trim
            ndims = length(size(adcp.(fld).(flds2{i})));
            trim_str = ['adcp.(fld).(flds2{i})(1:' num2str(n_trim)];
            for j = 1:(ndims-1)
                trim_str = [trim_str ',:'];
            end
            trim_str = [trim_str ')'];
            adcp.(fld).(flds2{i}) = eval(trim_str);
        end
    end
    adcp.(fld) = rmfield(adcp.(fld),'POS');
end



    