function ax = plot881aScan(fig,sonar)
% Function to quickly plot the two passes of a sonar 881a scan. This
% function relies on the sonar structure output from
% "parse881aDeployment.m" function.
%
% ax = plot881aScan(fig,sonar)
%
% Input
%   fig: figure on which to plot
%   sonar: single data structure containing time, angle, range, scan1,
%       scan2 fields
%
% Output
%   ax: axes object
%
% KJW
% 22 Jul 2026

% build cartesian grid
[R,PHI] = meshgrid(sonar.range,sonar.angle*pi/180);
PHI = -PHI; % flip because instrument is upside down on meltstake
X = R.*cos(PHI);
Y = R.*sin(PHI);

% set up axes
clf(fig)
pad = [.1 .1 .1 .1];
shift = [-.0 0];
for i = 2:-1:1
    ax(i) = axes(fig,'position',axgridpos(1,2,i,pad,shift),'color',[0 0 0]);
end

% plot
for i = 1:2
    pcolor(ax(i),X,Y,sonar.(['scan' num2str(i)]))
    shading(ax(i),'flat')
    axis(ax(i),'equal')
    clim(ax(i),extrema([extrema(sonar.scan1) extrema(sonar.scan2)]))
    xlabel(ax(i),'y [m]'); ylabel(ax(i),'x [m]');
    ax(i).Color = [0 0 0];
end
title(ax(1),sprintf('sonar 881a: %s',sonar.time))
colorbar(ax(2),'position',cbarpos(ax(2),.01,.02));