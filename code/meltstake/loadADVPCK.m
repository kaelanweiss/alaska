function ms = loadADVPCK(varargin)
% Function to load all ADV probe check output for all or selected meltstake
% deployments. The output is loaded into a structure that can be treated
% the same as the normal ADV structure in melt-related functions and
% scripts.
%
% ms = loadADVPCK()
% ms = loadADVPCK(deps)
% ms = loadADVPCK(tbl_path)
% ms = loadADVPCK(deps,tbl_path)
%
% Input
%   deps: (optional)
%   tbl_path: (optional)
%
% Output
%   ms: (1 x ndeps) structure containing ADV PCK output
%
% KJW
% 29 Jan 2024

% parse varargin
deps = nan;
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
switch nargin
    case 1
        if isstr(varargin{1})
            tbl_path = varargin{1};
        elseif isnumeric(varargin{1})
            deps = varargin{1};
        end
    case 2
        deps = varargin{1};
        tbl_path = varargin{2};
end

% read in table of deployment descriptions
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
mstbl = readtable(tbl_path,'sheet','meltcalculation');

% trim empty rows
mstbl = mstbl(~strcmp(mstbl.Folder,''),:);

% select deployments
if all(~isnan(deps))
    mstbl = mstbl(deps,:);
end

% grab deployment names and number of drill ops
depnames = mstbl.Folder;
ndrill = mstbl.Drill;
hann_width = mstbl.hann_width;
L = mstbl.L;
L0 = mstbl.L0;
x0 = [mstbl.x1 mstbl.x2 mstbl.x3];

% load
flds = {'pck_time','pck_dist','pck_amp','hdr_time','dp','dprobe_ave','dsvol_ave'};
ms(length(depnames)) = struct;

% loop through deployments
for i = 1:length(depnames)
    adv_path = fullfile('F:meltstake/data/raw',depnames{i},'adv','adv.mat');
    ms(i).depname = depnames{i};
    ms(i).ndrill = ndrill(i);
    ms(i).hann_width = hann_width(i);
    ms(i).L = L(i);
    ms(i).L0 = L0(i);
    ms(i).x0 = x0(i,:);
    % if file found, loop through fields
    if exist(adv_path,'file')
        load(adv_path,'adv')
        for j = 1:length(flds)
            ms(i).(flds{j}) = adv.(flds{j});
        end
    end
end

% % trim empty deployments
% no_pck = false(length(ms),1);
% for i = 1:length(ms)
%     if isempty(ms(i).pck_time)
%         no_pck(i) = true;
%     end
% end
% ms = ms(~no_pck);
