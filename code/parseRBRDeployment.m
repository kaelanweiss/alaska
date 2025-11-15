function rbr = parseRBRDeployment(depname,varargin)
% Function to collect RBR files and deployment information from a Dragon,
% Spider, etc. deployment into a single structure. Requires rks-tools
% library.
% 
% rbr = parseRBRDeployment(depname)
% rbr = parseRBRDeployment(depname,t1,t2)
% rbr = parseRBRDeployment(depname,rbr_dir)
% rbr = parseRBRDeployment(depname,t1,t2,rbr_dir)
%
% Input
%   depname: name of folder containing RBR/*.rsk files
%   t1: (optional) start of time window to parse
%   t2: (optional) end of time window to parse
%   rbr_dir: (optional) directory within depname in which to look for rsk
%   files
% 
% Output
%   rbr: 1 x (number of instruments) structure
%
% NOTE: 30 Jan 2024 changed datenum to datetime
%
% KJW
% Aug 2022

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
rbr_dir = 'rbr';
if nargin == 2 || nargin == 4
    rbr_dir = varargin{end};
end

% Get files
fnames = {dir(fullfile(depname,rbr_dir,'*.rsk')).name};
nf = length(fnames);

% Look for instrument.txt file
if exist(fullfile(depname,rbr_dir,'instruments.txt'),'file')
    fprintf('using instrument.txt file\n');
    rbr = parseInstrumentFile(depname,rbr_dir);
    sn = [rbr.sn];
    pos = [rbr.pos];
    
else
    fprintf('no instrument.txt file found\n');
    rbr(nf) = struct();
    sn = [];
end

% Check that directory exists
if ~exist(fullfile(depname,rbr_dir),'dir')
    warning('directory does not exist: %s\n',fullfile(depname,rbr_dir));
    return
end

% read files
for i = 1:nf
    fprintf('%d/%d\n',i,nf);

    % shift time vector
    if ~isfield(rbr(i),'t_offset') || isempty(rbr(i).t_offset) || isnan(rbr(i).t_offset)
        dt = 0;
    else
        dt = rbr(i).t_offset;
    end
%     if isfield(rbr(i),'t_offset') && (~isnan(rbr(i).t_offset) || ~isempty(rbr(i).t_offset))
%         dt = rbr(i).t_offset;
%     else
%         dt = seconds(0);
%     end

    rsk = RSKopen(fullfile(depname,rbr_dir,fnames{i}));
    rsk = RSKreaddata(rsk,'t1',datenum(t1)-dt/86400,'t2',datenum(t2)-dt/86400);
    % read data channels depending on what type of instrument
    if regexp(rsk.instruments.model,'solo')
        ch_idx = [1];
        modeli = 'solo';
    elseif regexp(rsk.instruments.model,'duet')
        ch_idx = [1 2];
        modeli = 'duet';
    elseif regexp(rsk.instruments.model,'concerto$')
        ch_idx = [1 2 3];
        modeli = 'concerto';
    elseif regexp(rsk.instruments.model,'concerto³')
        ch_idx = [1 2 3 5 6];
        modeli = 'concerto3';
    end
    rsk.data.values = rsk.data.values(:,ch_idx);
    ch_names = {rsk.channels(ch_idx).longName};
    ch_units = {rsk.channels(ch_idx).units};
    
    % find where to put the data using serial number
    sni = rsk.instruments.serialID;
    sn_idx = find(sni==sn);    
    if isempty(sn_idx) % if serial number not found in deployment, stick it at the end and fill in metadata
        if i == 1
            sn_idx = 1;
        else
            sn_idx = length(rbr)+1;
        end
        rbr(sn_idx).sn = sni;
        rbr(sn_idx).model = modeli;
        rbr(sn_idx).pos = NaN;
        rbr(sn_idx).t_offset = NaN;
    end
    
    % put data in output structure
    rbr(sn_idx).channels = ch_names;
    rbr(sn_idx).units = ch_units;
    rbr(sn_idx).time = datetime(rsk.data.tstamp,'convertfrom','datenum') + dt;
    rbr(sn_idx).values = rsk.data.values;
    
end

% sort by position if position is specified
if exist('pos','var')
    [~,sidx] = sort(pos);
    rbr = rbr(sidx);
end
    