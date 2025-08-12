% Script to plot profiles from a Dragon deployment.
%
% KJW
% 14 Sep 2022

clear;

% formatting
% set(0,'defaultTextInterpreter','latex');
% set(0,'defaultLegendInterpreter','latex');

% set berg name
bergs = {'20220822_oyster',...
         '20220823_queenofhearts',...
         '20220824_singingflower',...
         '20220825_singingflower'};
berg = bergs{4};

% set path to processed data
proc_path = 'F:/Alaska2022/data/iceberg_surveys/proc';
save_figs = false;
fig_path = fullfile(proc_path,berg,'dragon','figures');

% load "profs"
load(fullfile(proc_path,berg,'dragon','profiles.mat'));
n = length(profs);
pos = profs(1).rake.pos;

n = 10;
%% plot
% berg title
berg_splt = strsplit(berg,'_');
berg_ttl = sprintf('%s %s (local)',berg_splt{2},datestr(datenum(berg_splt{1},'yyyymmdd'),'dd-mmm'));

% general formatting
lw = 1.8;
fs = 16;

%%% just for cruise report %%%
tmp = profs(9).adcp.burst.vel;
bt = sigmaFilter(sigmaFilter(profs(9).adcp.bt.vel(2:end-1,:),1,1,1),2,1,1);
tmp = tmp - permute(repmat(bt,[1 1 85]),[1 3 2]);

beam_dirs = {'RD','RU','LU','LD'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 9%1:n
    figure(i); clf;
    ax = [];
    
    % rake temperature
    subplot(1,7,1:3); hold on;
    ax(end+1) = gca;
    pcolor([pos 2*pos(end)-pos(end-1)],profs(i).rake.depth,[profs(i).rake.temperature profs(i).rake.temperature(:,end)]);
    plot((profs(i).ctd.pos+0.5*diff(pos(1:2)))*[1 1],extrema(profs(i).rake.depth),'w--');
    shading flat;
    axis ij;
    ax1pos = get(gca,'position');
    cbar = colorbar;
    cbar.Position = [sum(ax1pos([1 3]))-.005 ax1pos(2) 0.01 ax1pos(4)];
    cbar.FontSize = fs-4;
    ylabel('depth (m)','fontsize',fs);
    xlabel('rake position (cm)','fontsize',fs)
    title([berg_ttl ' |  Dragon profile ' num2str(i)],'fontsize',fs+2);
    cmocean('thermal');
    set(gca,'XTick',pos+0.5*diff(pos(1:2)));
    set(gca,'XTickLabel',strsplit(num2str(pos),' '));
    xlim(extrema([profs(i).rake.pos(1) 2*pos(end)-pos(end-1)]))
    ylim(extrema(profs(i).rake.depth))
    text(pos(end-1),2,'T ($^{o}$C)','fontsize',fs)
    
    % adcp
    p0 = 0.3;
    for j = 1:4
        subplot(1,7,j+3);
        ax(end+1) = gca;
        pcolor(profs(i).adcp.burst.range,profs(i).adcp.burst.pressure-p0,tmp(:,:,j))
        shading flat
        axis ij
        cmocean('balance')
        set(gca,'yaxislocation','right')
        set(gca,'xdir','reverse')
        caxis(0.1*[-1 1])
        xlabel('beam dist (m)','fontsize',fs)
        title(['beam: ',beam_dirs{j}],'fontsize',fs)
        if j == 4
            axpos = get(gca,'position');
            set(gca,'ytick',[])
            cbar = colorbar;
            cbar.Position = [axpos(1)+axpos(3) axpos(2) .01 axpos(4)];
            cbar.Label.String = 'beam velocity (m/s)';
            cbar.Label.FontSize = fs-4;
            cbar.FontSize = fs-4;
        end
    end    
    linkaxes(ax,'y');
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
    
    