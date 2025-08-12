% Script to create section comparison plots to dive into interesting
% section pairs. These are figures 4 and 5 of the masters thesis.
%
% KJW
% 8 Apr 2024
set(0,'defaulttextinterpreter','latex');

clear

% load data
proc_path = 'F:meltstake\data\proc';
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

% define pairs by deployment index (table row)
pairs = [9 15;...
         8 29;...
         9 22;...
         9 21;...
         9 3;...
         9 18];%...
%          25 31;...
%          26 31];

% zoom into time limits
t0 = [datetime(2023,7,9,1,40,0) datetime(2023,7,10,19,30,0);...
      datetime(2023,7,8,22,50,0) datetime(2023,9,23,21,38,0);...
      datetime(2023,7,9,1,40,0) datetime(2023,7,11,1,4,0);...
      datetime(2023,7,9,1,40,0) datetime(2023,7,11,1,8,0);...
      datetime(2023,7,9,1,22,0) datetime(2023,5,29,21,18,0);...
      NaT NaT];
t_width = [minutes(15);...
           minutes(2);...
           minutes(26);...
           minutes(30);...
           minutes(10)];

% melt rate
m = msTable2Vector(ms_tbl.m);
m_ci = ms_tbl.m_ci;
m = m(pairs)*.24;
m_ci = m_ci(pairs)*.24;

% preallocate data space
np = size(pairs,1);
T_all = cell(np,2);
adcp_all = cell(np,2);
uT_all = cell(np,2);
m_all = nan(np,2);
dep_nums = nan(np,2);
wind_nums = nan(np,2);

% hanning width
hann_dt = 3; % s

for i = 1:np    
    % load section data
    for j = 1:2
        dep_num = ms_tbl.Number(pairs(i,j));
        wind_num = ms_tbl.Window(pairs(i,j));
        fprintf('%d.%d\n',dep_num,wind_num)
        
        load(fullfile(proc_path,ms_tbl.Folder{pairs(i,j)},sprintf('T%d.mat',wind_num)))
        load(fullfile(proc_path,ms_tbl.Folder{pairs(i,j)},sprintf('adcp%d.mat',wind_num)))

        % flip order of thermistors if need be
        if any(pairs(i,j)==[15])
            T.T = fliplr(T.T);
        end
        
        % convert ADCP range to distance from ice (estimate)
        rmax = ms_tbl.rmax(pairs(i,j));
        adcp.burst.range = rmax/cosd(25)-adcp.burst.range;% + 0.5*adcp.burst.cellsize;

        % calculate velocity scale
        vel = adcp.burst.vel_ice(:,:,1:2); % u,w
        % smooth across 3 points in the profile
        for k = 1:size(vel,1)
            for p = 1:size(vel,3)
                vel(k,:,p) = hannFilter(squeeze(vel(k,:,p)),3);
            end
        end
        % calculate magnitude
        vel_mag = vecnorm(vel,2,3);
        % profile max and mean
        vel_max = max(vel_mag,[],2);
        adcp.burst.vel_scale = vel_max;

        % calculate velocity profiles
        adcp.burst.vel_prof = squeeze(mean(cat(3,vel,vel_mag),1,'omitnan'));
        adcp.burst.vel_prof_std = squeeze(std(cat(3,vel,vel_mag),1,'omitnan'));

        % u*T
        vel_interp = interp1(adcp.burst.time,adcp.burst.vel_scale,T.time,'nearest',nan);
        T_mean = mean(T.T,2);
        uT = (vel_interp/mean(vel_interp,'omitnan')).*(T_mean/mean(T_mean,'omitnan'));

        % smoothing
        k_T = round(hann_dt/seconds(diff(T.time(1:2))));
        k_u = hann_dt*adcp.burst.samplerate;
        % T
        for k = 1:size(T.T,2)
            T.T(:,k) = hannFilter(T.T(:,k),k_T);
        end
        % U
        adcp.burst.vel_scale = hannFilter(adcp.burst.vel_scale,k_u);
        for k = 1:size(adcp.burst.vel,2)
            for p = 1:2
                adcp.burst.vel_ice(:,k,p) = hannFilter(adcp.burst.vel_ice(:,k,p),k_u);
            end
        end

        T_all{i,j} = T;
        adcp_all{i,j} = adcp;
        uT_all{i,j} = uT;
        dep_nums(i,j) = dep_num;
        wind_nums(i,j) = wind_num;
    end
end

%% plot
% plot params
% padding
pxi = .05;
pxo = .12;
pyi = .02;
pyo = .08;

lw = 1.2;
fs = 11;

% T colors
T_clr = [0.2 1 0.5]/1.2;
T_clrs = T_clr'*[1 0.5 0];
T_lgd_lbls = {'near','mid','far'};

% labels
panel_lbls = {'$T$','max($|u|$)','$u$','$w$'};
y_lbls = {'$^\circ$C','m/s',{'dist from','ice (m)'},{'dist from','ice (m)'}};

% pcolor lims
pc_lims = [.1 .1 .2 .2 .1 .1];

for i = 1:np
    % figure i
    % create and size figure
    figsize = [18 12];
    fig = figure(i); clf
    set(fig,'units','centimeters'); 
    fpos = get(fig,'position');
    fpos(3:4) = figsize;
    set(fig,'position',fpos);
    set(fig,'paperunits','centimeters')
    set(fig,'papersize',figsize)
    clear ax
    
    % axis limits
    vel_lim = max([extrema(adcp_all{i,1}.burst.vel_scale) extrema(adcp_all{i,2}.burst.vel_scale)]);
    T_lim = extrema([extrema(T_all{i,1}.T) extrema(T_all{i,2}.T)]) + 0.5*[-1 1];
    u_lim = max(abs(extrema([extrema(adcp_all{i,1}.burst.vel_ice(:,:,1)) extrema(adcp_all{i,2}.burst.vel_ice(:,:,1))])));
    w_lim = max(abs(extrema([extrema(adcp_all{i,1}.burst.vel_ice(:,:,2)) extrema(adcp_all{i,2}.burst.vel_ice(:,:,2))])));
    uw_lim = max([u_lim w_lim]);

    % loop through each section
    for j = 1:2
        % plot
        % T
        ax(4*(j-1)+1) = axes(fig,'position',axgridpos(4,2,4*(j-1)+1,pxi,pyi,pxo,pyo,'flip'));
        hold on
        for k = 1:3
            plot(T_all{i,j}.time,T_all{i,j}.T(:,k),'linewidth',lw,'color',T_clrs(:,k))
        end
        title(sprintf('%s (m: %.1f$\\pm$%.1f m/day)',char(64+2*(i-1)+j),m(i,j),m_ci(i,j)),'fontsize',fs)
        ylim(T_lim)
        box on

        % vel scale
        ax(4*(j-1)+2) = axes(fig,'position',axgridpos(4,2,4*(j-1)+2,pxi,pyi,pxo,pyo,'flip'));
        plot(adcp_all{i,j}.burst.time,hannFilter(adcp_all{i,j}.burst.vel_scale,8),'k','linewidth',lw)
        ylim([0 vel_lim])

        % u
        ax(4*(j-1)+3) = axes(fig,'position',axgridpos(4,2,4*(j-1)+3,pxi,pyi,pxo,pyo,'flip'));
        pcolor(adcp_all{i,j}.burst.time,adcp_all{i,j}.burst.range,adcp_all{i,j}.burst.vel_ice(:,:,1)')
        shading flat
        cmocean('bal')
        clim(pc_lims(i)*[-1 1])

        % w
        ax(4*(j-1)+4) = axes(fig,'position',axgridpos(4,2,4*(j-1)+4,pxi,pyi,pxo,pyo,'flip'));
        pcolor(adcp_all{i,j}.burst.time,adcp_all{i,j}.burst.range,adcp_all{i,j}.burst.vel_ice(:,:,2)')
        shading flat
        cmocean('bal')
        clim(pc_lims(i)*[-1 1])

        % link axes
        linkaxes(ax((1:4)+4*(j-1)),'x')
        xlim(ax(4*j),extrema(T_all{i,j}.time))

    end

    % clear internal x ticks
    set(ax(setxor(1:8,[4 8])),'xticklabels',{})

    % colorbars
    cbu = colorbar(ax(7));
    cbw = colorbar(ax(8));
    cbu.Position = cbarpos(ax(7),.005,.01);
    cbw.Position = cbarpos(ax(8),.005,.01);
    cbw.Position(4) = sum(cbu.Position([2 4])) - cbw.Position(2);
%     cbu.Label.String = 'm/s';
%     cbu.Label.Interpreter = 'latex';
%     cbu.Label.FontSize = fs;
    delete(cbu);
    cbw.Label.String = 'm/s';
    cbw.Label.Interpreter = 'latex';
    cbw.Label.FontSize = fs;
    
    % panel and y axis labels
    for j = 1:4
        for k = 1:2
            txt = text(ax(j+4*(k-1)),.02,.98,['(' char(96 + 4*(k-1) + j) ') ' panel_lbls{j}],'units','normalized','fontsize',fs,...
            'verticalalignment','top','horizontalalignment','left');
            if j >= 3
%                 txt.BackgroundColor = 'w';
            end
        end
        ylabel(ax(j),y_lbls{j},'fontsize',fs)
    end 
    
    % adjust time windows
    if ~isnat(t0(i,1))
        xlim(ax(1),t0(i,1)+[0 t_width(i)])
        xlim(ax(5),t0(i,2)+[0 t_width(i)])
    end

    % adjust ADCP pcolor distance ranges
    for j = [3 4 7 8]
        ylim(ax(j),[0 max([adcp_all{i,1}.burst.range(1) adcp_all{i,2}.burst.range(1)])])
    end

    % T legend
    legend(ax(1),T_lgd_lbls,'orientation','horizontal','location','south',...
        'fontsize',fs-1,'edgecolor','w')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % figure i + 10
%     fig = figure(i+10); clf
%     clear ax
% 
%     % loop through each section
%     for j = 1:2
%         % plot
%         % vel profs
%         ax(2*(j-1)+1) = axes(fig,'position',axgridpos(2,2,2*(j-1)+1,pxi,4*pyi,pxo,1.2*pyo,'flip'));
%         hold on
%         title(sprintf('%d.%d (%.2f$\\pm$%.2f m/day)',dep_nums(i,j),wind_nums(i,j),m(i,j),m_ci(i,j)))
%         plot([0 2],[0 0],'k:','linewidth',lw)
%         clrs = {colors(1),colors(2),'k'};
%         for k = 1:3
%             errorbar(adcp_all{i,j}.burst.range,adcp_all{i,j}.burst.vel_prof(:,k),adcp_all{i,j}.burst.vel_prof_std(:,k),...
%                 'color',clrs{k},'linewidth',lw,'marker','.')
%         end
%         xlim([0 max([adcp_all{i,1}.burst.range(1) adcp_all{i,2}.burst.range(1)])])
%         xlabel('adcp range (m)','fontsize',fs)
% 
%         % u*T histogram
%         ax(2*(j-1)+2) = axes(fig,'position',axgridpos(2,2,2*(j-1)+2,pxi,4*pyi,pxo,1.2*pyo,'flip'));
%         histogram(uT_all{i,j},25)
%         xlabel('$u$*$T$/$\overline{u}\overline{T}$')
% 
%     end
% 
%     % match plot limits between pairs
%     linkaxes(ax([1 3]),'y')
%     linkaxes(ax([2 4]),'x')
%     
%     % panel and y axis labels
%     ylabel(ax(1),'vel profs (m/s)','fontsize',fs)

end
