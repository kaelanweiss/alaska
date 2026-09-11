% script to make gif from sonar 881a deployment
%
% KJW
% 6 Aug 2026

clear

raw_dir = 'F:\GIANT2026\meltstake\data\raw';
fig_dir = 'F:\GIANT2026\meltstake\figures';
depname = '20260730_1530_ms01';
% depname = '20260803_1430_ms01';
% depname = '20260816_1850_ms01';
load(fullfile(raw_dir,depname,'sonar881a\sonar881a.mat'))

%% set up figure
ns = length(sonar);
% build cartesian grid
[R,PHI] = meshgrid(sonar(1).range,sonar(1).angle*pi/180);
PHI = -PHI; % flip because instrument is upside down on meltstake
X = R.*cos(PHI);
Y = R.*sin(PHI);


fig = figure(1); clf
setFigureSize(fig,[14 9])
ax = axes(fig);
ax.Position(3) = ax.Position(3)*0.95;

h = pcolor(Y',X',sonar(1).scan1');
shading flat
cmocean('-mat')
ax.Color = 'k';
axis(ax,'equal')
axpos = ax.Position;
cbar = colorbar(ax,'position',cbarpos(ax,.01,.025));
ax.Position = axpos;
cbar.Label.String = 'amplitude';
cbar.Label.FontSize = 11;
xlabel('along-ice (m)')
ylabel('ice-normal (m)')

ttl_str = {[strrep(depname,'_','\_') ' | '],'timestamp'};
ttl = title([ttl_str{1} ttl_str{2}]);

%%
amp = nan([ns size(sonar(1).scan1)]);
for i = 1:ns
    amp(i,:,:) = flipud((sonar(i).scan1+sonar(i).scan2)/2);
end

% for j = 1:size(amp,2)
%     fprintf('%d/%d\n',j,size(amp,2))
%     for k = 1:size(amp,3)
%         amp(:,j,k) = hannFilter(amp(:,j,k),5);
%     end
% end

%% create images
im = cell(1,ns);
fprintf('Creating images\n')
for i = 1:ns
i = 1;
    fprintf('%d/%d\n',i,ns)
    h.CData = squeeze(amp(i,:,:))';
    ttl.String = [ttl_str{1} sprintf('%s',sonar(i).time)];
    ax.CLim = [0 80];
    drawnow
    frame = getframe(fig);
    im{i} = frame2im(frame);
end

%% write file
skip = 1;
dt = 30/(ns/skip);
if ~exist(fullfile(fig_dir,depname),'dir')
    mkdir(fullfile(fig_dir,depname))
end

filename = fullfile(fig_dir,depname,['sonar881a_' depname '.gif']);
fprintf('Writing frames to file\n')
for i = 1:skip:ns
    fprintf('%d/%d\n',i,ns)
    [A,map] = rgb2ind(im{i},256);
    if i == 1
        imwrite(A,map,filename,'gif','loopcount',Inf,'delaytime',dt)
    else
        imwrite(A,map,filename,'gif','writemode','append','delaytime',dt)
    end
end