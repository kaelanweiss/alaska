function [ax,cbar] = adcpQuickPlot(fig,adcp,fld,cbarlim,tlim,rlim,dt,varargin)
% Function to plot adcp fields as panel of pcolors. This is meant to
% streamline the process of making tons of multi-panel figures for an ADCP
% deployment.
%
% [ax,cbar] = adcpQuickPlot(fig,adcp,fld,cbarlim,tlim,rlim,dt)
% [ax,cbar] = adcpQuickPlot(fig,adcp,fld,cbarlim,tlim,rlim,dt,subfld)
%
% Input
%   fig: figure handle (e.g. figure(1))
%   adcp: adcp data structure (e.g. output from parseNortekDeployment)
%   fld: field to plot within the "burst" substructure (e.g. "vel", "amp")
%   cbarlim: limits of colorbar as vector [x y]
%   tlim: time limits for plotting as vector of datetimes [t1 t2], input
%       NaT to plot all data
%   rlim: range limit as vector, input NaN to plot all data
%   dt: time steps to skip if plotting lots of data, i.e. 10 plots every 10
%       points in time
%   subfld: (optional) choose a subfield within the adcp structure other
%       than the "burst" subfield e.g. "echo"
%
% KJW
% 2023

if (isdatetime(tlim) && all(isnat(tlim))) || (~isdatetime(tlim) && all(isnan(tlim)))
    tlim = [datetime([0 1 1]) datetime([3000 1 1])];
end

if isnan(rlim)
    rlim = [-5 Inf];
end

switch fld(1:3)
    case 'vel'
        cmap = 'bal';
        units = 'm/s';
    case 'cor'
        cmap = 'mat';
        units = '%';
    case 'amp'
        cmap = 'thermal';
        units = 'dB';
end

subfld = 'burst';
if nargin > 7
    subfld = varargin{1};
end

[~,~,nb] = size(adcp.(subfld).(fld));
t = adcp.(subfld).time;
if ~isdatetime(t)
    t = datetime(t,'convertfrom','datenum');
end

idxt = find(t>=tlim(1) & t<tlim(2));
idxt = idxt(1:dt:end);
idxr = adcp.(subfld).range>=rlim(1) & adcp.(subfld).range<=rlim(2);

t = t(idxt);
r = adcp.(subfld).range(idxr);

fs = 12;
fig;
clf(fig)
clear ax
for i = 1:nb
    axpos = axgridpos(nb,1,i,0,0.02,0.085,0.07);
    ax(i) = axes(fig,'position',axpos);
    pcolor(ax(i),t,r,adcp.(subfld).(fld)(idxt,idxr,i)')
    shading(ax(i),'flat')
    clim(ax(i),cbarlim)
    cmocean(cmap,ax(i))
    ax(i).FontSize = 9;
    ylabel('range (m)','fontsize',fs)
    if max(r) > 3
        axis ij
    end
    if i<nb
        set(ax(i),'xticklabel','')
    else
        datetick(ax(i),'x','mmmdd HH:MM','keeplimits')
    end
    box on
end
cbar = colorbar(ax(nb));
ax(nb).Position = axpos;
cbar.Position = [sum(ax(nb).Position([1 3]))+0.005 ax(nb).Position(2) 0.008 sum(ax(1).Position([2 4]))-ax(nb).Position(2)];
cbar.Label.String = sprintf('%s (%s)',strrep(fld,'_','\_'),units);
cbar.Label.FontSize = fs;

linkaxes(ax)
fig.Color = 0.85*[1 1 1];