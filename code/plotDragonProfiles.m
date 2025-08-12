% Script to plot profiles from a Dragon deployment.
%
% KJW
% 14 Sep 2022

%clear;

% formatting
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

% set berg name
bergs = {'20220820_marchhare',...
         '20220821_marchhare',...
         '20220822_oyster',...
         '20220823_queenofhearts',...
         '20220824_singingflower',...
         '20220825_singingflower'};
berg = bergs{5};

% set path to processed data
proc_path = 'F:/Alaska2022/data/iceberg_surveys/proc';
save_figs = false;
fig_path = fullfile(proc_path,berg,'dragon','figures');

% load "profs"
load(fullfile(proc_path,berg,'dragon','profiles.mat'));
n = length(profs);
pos = profs(1).rake.pos;

%% plot
% berg title
berg_splt = strsplit(berg,'_');
berg_ttl = sprintf('%s %s (local)',berg_splt{2},datestr(datenum(berg_splt{1},'yyyymmdd'),'dd-mmm'));
ctd_depth_offset = 0.35;

% general formatting
lw = 1.5;
fs = 14;
for i = 1:n
    figure(i); clf;
    clear ax

    fig_ttl = {[berg_ttl ' |  profile ' num2str(i)], [' ' datestr(profs(i).ctd.time(1),'(mmmdd HH:MM:SSZ)')]};
    
    % rake temperature
    subplot(1,7,1:3); hold on;
    ax(1) = gca;
    ax(1).Position(3) = ax(1).Position(3)-0.015;
    pcolor([pos 2*pos(end)-pos(end-1)],profs(i).rake.depth,[profs(i).rake.temperature profs(i).rake.temperature(:,end)]);
    plot((profs(i).ctd.pos+0.5*diff(pos(1:2)))*[1 1],extrema(profs(i).rake.depth),'w--');
    shading flat;
    axis ij;
    ax1pos = get(gca,'position');
    cbar = colorbar;
    cbar.Position = [sum(ax1pos([1 3]))+.0025 ax1pos(2) 0.01 ax1pos(4)];
    ylabel('depth (m)','fontsize',fs);
    xlabel('rake position (cm)','fontsize',fs)
    %title([berg_ttl ' |  profile ' num2str(i)],'fontsize',fs+2);
    text(0,1.01,fig_ttl,'fontsize',fs+2,'units','normalized','VerticalAlignment','baseline')
    text(1.01,1,'T (\circC)','fontsize',fs,'units','normalized',...
        'VerticalAlignment','bottom','HorizontalAlignment','center','Interpreter','tex')
    cmocean('thermal');
    set(gca,'XTick',pos+0.5*diff(pos(1:2)));
    set(gca,'XTickLabel',strsplit(num2str(pos),' '));
    xlim(extrema([profs(i).rake.pos(1) 2*pos(end)-pos(end-1)]))
    ylim(extrema(profs(i).rake.depth))
    if diff(clim) > 2
        clim(clim + 0.25*[1 -1])
    end
    clim([5.5 7])
    
%     % ctd temperature
%     ctd_depth_offset = 0.2;
%     subplot(1,7,4);
%     ax(end+1) = gca;
%     plot(profs(i).ctd.temperature,profs(i).ctd.depth-ctd_depth_offset,'-','linewidth',lw,'color',colors(1));
%     axis ij;
%     grid;
%     box on;
%     set(gca,'YTickLabel',{});
%     xlabel('CTD temp (C)','fontsize',fs);
%     xlim(extrema(profs(i).ctd.temperature)+[-.2 .2]);
%     
%     % ctd salinity
%     subplot(1,7,5);
%     ax(end+1) = gca;
%     plot(profs(i).ctd.salinity,profs(i).ctd.depth-ctd_depth_offset,'-','linewidth',lw,'color',colors(2));
%     axis ij;
%     grid;
%     box on;
%     set(gca,'YTickLabel',{});
%     xlabel('CTD sal (psu)','fontsize',fs);
%     xlim(extrema(profs(i).ctd.salinity)+[-.2 .2]);

    % ctd T and S
    subplot(1,7,4:5)
    ax(2) = gca;
    % T
    plot(profs(i).ctd.temperature,profs(i).ctd.depth-ctd_depth_offset,'-','linewidth',lw,'color',colors(1));
    axis ij;
    grid;
    box on;
    set(gca,'YTickLabel',{});
    xlabel('CTD temp (C)','fontsize',fs);
    xlim(extrema(profs(i).ctd.temperature)+[-.2 .2]);
    ax(2).XColor = colors(1)*0.7;
    box off

    % S
    ax(3) = axes(figure(i),'position',ax(2).Position);
    plot(profs(i).ctd.salinity,profs(i).ctd.depth-ctd_depth_offset,'-','linewidth',lw,'color',colors(2));
    axis ij;
    grid;
    box on;
    set(gca,'YTickLabel',{});
    xlabel('CTD sal (psu)','fontsize',fs);
    xlim(extrema(profs(i).ctd.salinity)+[-.2 .2]);
    ax(3).XAxisLocation = 'top';
    ax(3).Color = 'none';
    ax(3).XColor = colors(2)*0.85;
    box off

    % heading
    subplot(1,7,6);
    ax(4) = gca;
    plot(profs(i).rov.heading,profs(i).rov.depth,'k-','linewidth',lw);
    axis ij;
    grid;
    set(gca,'YTickLabel',{});
    xlabel('heading (deg)','fontsize',fs);
    xlim(extrema(profs(i).rov.heading));
    
    % roll/pitch
    subplot(1,7,7);
    ax(5) = gca;
    hold on;
    plot(profs(i).rov.pitch,profs(i).rov.depth,'k-','linewidth',lw);
    plot(profs(i).rov.roll,profs(i).rov.depth,'k:','linewidth',lw);
    axis ij;
    grid;
    box on;
    set(gca,'YAxisLocation','right')
    xlabel('pitch/roll (deg)','fontsize',fs);
    legend({'pitch','roll'},'fontsize',fs-2,'location','southwest');
    ylabel('depth (m)','fontsize',fs);
    
    linkaxes(ax,'y');
    ylim(extrema(profs(i).rov.depth));

end

%%
% dur = 8;
% for i = 1:n
%     figure(i);
%     fprintf('%d) ',i);
%     for j = 1:dur
%         fprintf('%d...',dur+1-j)
%         pause(1)
%     end
%     fprintf('0\n')
% end

%%
if save_figs
    if ~exist(fig_path,'dir')
        mkdir(fig_path);
    end
    
    for i = 1:n
        fprintf('saving figure %d\n',i);
        print(figure(i),fullfile(fig_path,[berg '_profile' sprintf('%02d',i) '.png']),'-dpng','-r400');
    end
end
    
    