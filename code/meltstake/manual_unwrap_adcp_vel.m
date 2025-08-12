% Script to manually unwrap ADCP velocities by inspection. This has some
% obvious limitations, but may work in some simple cases.
%
% KJW
% 28 Aug 2024

clear

depname = 'ms03_20230709_2138';

raw_path = fullfile('F:meltstake/data/raw/',depname,'adcp');
load(fullfile(raw_path,'adcp.mat'));

%%
% max range
rmax = 0.6;
idxr = adcp.burst.range<=rmax;
range = adcp.burst.range(idxr);

% data dimensions
[nt,~,nb] = size(adcp.burst.vel);
nc = sum(idxr);

% vel range
vmax = max(abs(extrema(adcp.burst.vel)));
v_amb = diff(extrema(adcp.burst.vel(:,:,1:4)))*ones(1,nb);
v_amb(5) = diff(extrema(adcp.burst.vel(:,:,5)));

% plot data
axv = adcpQuickPlot(figure(1),adcp,'vel',vmax,NaT,[0 rmax],1);
axc = adcpQuickPlot(figure(3),adcp,'cor',[40 100],NaT,[0 rmax],1);

% QC
cor_min = 50;
vel = adcp.burst.vel;
vel(adcp.burst.cor<cor_min) = nan;

% wrapped intervals
wrapped = cell(nb,nc);
if ~exist(fullfile(raw_path,'wrapped.mat'),'file')
    temp = [];
    save(fullfile(raw_path,'wrapped.mat'),'temp')
end

% clicking figure
figure(2); clf
clear ax2
ax2(1) = subplot(2,1,1);
cmocean('bal')
ax2(2) = subplot(2,1,2);
linkaxes(ax2,'x')

%% find wrapped intervals manually
for i = 1%:nb
    % plot upper panel
    cla(ax2(1)); hold(ax2(1),'on')
    pcolor(ax2(1),1:nt,range,vel(:,idxr,i)')
    shading(ax2(1),'flat')
    clim(vmax*[-1 1])
    pbox = plot(ax2(1),1:nt,[range(1) range(1)+adcp.burst.cellsize]'*ones(1,nt),'k--');

    for j = 1%:nc
        fprintf('---b%dc%02d---\n',i,j)
        % update upper panel lines
        pbox(1).YData = adcp.burst.range(j)*ones(1,nt);
        pbox(2).YData = adcp.burst.range(j+1)*ones(1,nt);
        
        % plot lower panel
        cla(ax2(2)); hold(ax2(2),'on')
        plot(ax2(2),1:nt,vel(:,j,i),'.-')
        plot(ax2(2),1:nt,0.5*v_amb(i)*[-1 1]'*ones(1,nt),'k--')
        xlim([1 nt])

        % setup
        intervals = [];

        s = '';
        while ~strcmp(s,'stop')
            s = input('','s');
            [idx_get,~] = ginput(2);
            if length(idx_get)==2
                intervals = cat(1,intervals,idx_get');
            else
                warning('single point')
            end
        end
        
        % save
%         wrapped{i,j} = intervals;
%         varname = sprintf('b%dc%02d',i,j);
%         eval(sprintf('%s = intervals;',varname))
%         save(fullfile(raw_path,'wrapped.mat'),varname,'-append')
    end
end

%% fix the velocities
vel_fix = vel;

% beams
for i = 1:nb
    % cells
    for j = 1:nc
        intervals = wrapped{i,j};
        n_int = size(intervals,1);
        
        % wrapped intervals
        for k = 1:n_int
            idx1 = ceil(intervals(k,1));
            idx2 = floor(intervals(k,2));
            
            % check that the interval isn't a double-click
            if idx2-idx1 >= 0
                v_slice = vel(idx1:idx2,j,i);
                % push up or down
                if mean(v_slice,'omitnan')<0
                    vel_fix(idx1:idx2,j,i) = v_slice + v_amb(i);
                elseif mean(v_slice,'omitnan')>0
                    vel_fix(idx1:idx2,j,i) = v_slice - v_amb(i);
                elseif mean(v_slice,'omitnan')==0
                    warning('v_slice mean = 0 (b%dc%02d)',i,j)
                end
            end
        end
    end
end

adcp.burst.vel_unwrap = vel_fix;

axf = adcpQuickPlot(figure(4),adcp,'vel_unwrap',2*vmax,NaT,NaN,1);

%% inspect fixes
figure(5); clf
hold on
pc1 = plot(1:nt,zeros(1,nt),'.-','color',colors(2));
pc2 = plot(1:nt,zeros(1,nt),'.-','color',colors(1));
plot(1:nt,0.5*v_amb(1)*[-1 1]'*ones(1,nt),'k--')

for i = 1:nb
    for j = 1:nb
        title(sprintf('beam %d cell %d',i,j))
        pc1.YData = vel_fix(:,j,i);
        pc2.YData = vel(:,j,i);
        input(' ')
    end
end

%% inspect fixes v2
figure(6); clf
clear ax6
for i = 1:nb
    ax6(1) = subplot(2,1,1);
    pcolor(1:nt,range,vel(:,idxr,i)')
    shading flat
    cmocean('bal')
    colorbar
    clim(0.2*[-1 1])

    ax6(2) = subplot(2,1,2);
    pcolor(1:nt,range,vel_fix(:,idxr,i)')
    shading flat
    cmocean('bal')
    colorbar
    clim(0.2*[-1 1])

    linkaxes(ax6)
    ylim([0.12 0.63])

    title(ax6(1),sprintf('beam %d',i))

    input('')
end



