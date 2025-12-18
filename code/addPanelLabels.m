function lbls = addPanelLabels(ax,fs,varargin)
% Function to add text labels in upper left corner of axes.
%
% addPanelLabels(ax,fs)
% addPanelLabels(ax,fs,fmt)
% addPanelLabels(ax,fs,[],flip,nrows)
%
% Input
%   ax: list of axes
%   fs: fontsize
%   fmt: (optional) format string, must include %c, default is '%c)'
%   flip: (optional) input 'flip' to count column-wise
%   nrows: (optional) number of rows, must be included if 'flip' is
%          specified
%   
% Output
%   lbls: text objects
%
% KJW
% 17 Dec 2025

% parse input
p = inputParser;
validStrInput = @(x) isstring(x) || ischar(x);
validFmtInput = @(x) isempty(x) || (validStrInput(x) && length(regexp(x,'%c')) == 1);
addRequired(p,'ax');
addRequired(p,'fs');
addOptional(p,'fmt','%c)',validFmtInput);
addOptional(p,'flip','noflip',validStrInput);
addOptional(p,'nrows',1);
parse(p,ax,fs,varargin{:});
% fmt
fmt = p.Results.fmt;
if isempty(fmt)
    fmt = '%c)';
end
% flip
flip = p.Results.flip;
if strcmp(flip,'noflip')
    flip = 0;
elseif ~strcmp(flip,'noflip') && strcmp(p.Results.flip,'flip')
    flip = 1;
    nrows = p.Results.nrows;
end

% flip the axes
n = length(ax);
idx = 1:n;
if flip
    grid = reshape(idx,[nrows n/nrows])';
    idx = grid(:)';
end
for i = 1:n
    lbls(i) = text(ax(idx(i)),.02,1,sprintf(fmt,96+i),'units','normalized','fontsize',fs,'horizontalalignment','left','verticalalignment','top');
end
a=1;


