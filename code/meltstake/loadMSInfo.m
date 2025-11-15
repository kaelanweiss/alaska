function mstbl = loadMSInfo(varargin)
% Function to load meltstake deployment information from the spreadsheet
% "meltstake_deployments.xlsx". This can be useful for loading and working
% with all or a subset of deployments. Different sheets within the
% speadsheet can be specified. This is basically just a wrapper for the
% "readtable" function.
%
% ms = loadMSInfo()
% ms = loadMSInfo(deps)
% ms = loadMSInfo(sheet)
% ms = loadMSInfo(deps,sheet)
%
% Input
%   deps: (optional) deployment number(s) as they appear in the "Number"
%       column, default is all deployments
%   sheet: (optional) name of sheet to load, default is 'overview'
%
% Output
%   ms: table containing deployment information
%
% KJW
% 30 Jan 2024

TBL_PATH = 'G:\Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';

% parse varargin
deps = nan;
sheet = 'overview';
switch nargin
    case 1
        if isstr(varargin{1})
            sheet = varargin{1};
        elseif isnumeric(varargin{1})
            deps = varargin{1};
        end
    case 2
        deps = varargin{1};
        sheet = varargin{2};
end

% read in table of deployment descriptions
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
mstbl = readtable(TBL_PATH,'sheet',sheet);

% trim empty rows
mstbl = mstbl(~strcmp(mstbl.Folder,''),:);

% select deployments
if isnan(deps)
    return
else
    dep_idx = false(size(mstbl,1),1);
    for i = 1:length(deps)
        dep_idx = dep_idx | mstbl.Number == deps(i);
    end
    mstbl = mstbl(dep_idx,:);
end

