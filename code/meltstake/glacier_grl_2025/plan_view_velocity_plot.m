% Script to plot plan view distributions of ocean velocity from Polly/Aries 
% to contextualize MS observations.
%
% KJW
% 10 Dec 2025
clear

% load
load platform_comparison_data\grid_vel.mat

%% plot
fs = 11;
q_scale = 4;
lw = 1;
fig_size = [14 8];
fig = figure(1); clf
setFigureSize(fig,fig_size);

leg_vec = 0.5*[-1 1]/sqrt(2);
leg_vec_pos = [56.835 -132.359];

for i = 1:3
    dlat = mean(diff(grid_vel(i).lat));
    dlon = dlat/cosd(mean(grid_vel(i).lat));
    
    % plot
    ax = subplot(1,3,i);
    hold on
    box on
    
    % terminus line
    h(1) = patch(grid_vel(i).term.lon,grid_vel(i).term.lat,0.2*[1 1 1],'facealpha',0.25,'edgecolor','k','linewidth',lw);
      
    % velocity arrows
    h(4) = quiver(grid_vel(i).lon(1:end-1),grid_vel(i).lat(1:end-1),grid_vel(i).vel(:,:,1,2)*q_scale*dlon,grid_vel(i).vel(:,:,2,2)*q_scale*dlat,'off','r','linewidth',lw);
    h(3) = quiver(grid_vel(i).lon(1:end-1),grid_vel(i).lat(1:end-1),grid_vel(i).vel(:,:,1,1)*q_scale*dlon,grid_vel(i).vel(:,:,2,1)*q_scale*dlat,'off','k','linewidth',lw);
    
    % ms location
    h(2) = plot(grid_vel(i).ms_loc(2),grid_vel(i).ms_loc(1),'ro','markersize',4,'markerfacecolor','r','linewidth',lw);
    
    legend(h,{'glacier','ms location','sfc',sprintf('ms depth (%d m)',round(grid_vel(i).depth{1}(2)))},'location','northeast','autoupdate','off','fontsize',fs-3)
    quiver([0 leg_vec_pos(2)],[0 leg_vec_pos(1)],[0.05 1]*leg_vec(1)*q_scale*dlon,[0.05 1]*leg_vec(2)*q_scale*dlat,'off','k','linewidth',lw);
    text(leg_vec_pos(2),leg_vec_pos(1)+leg_vec(2)*q_scale*dlat,sprintf('%.1f m/s',round(vecnorm(leg_vec),1)),'fontsize',fs-2,'horizontalalignment','center','verticalalignment','bottom')
    xlim(-132-[0.367 0.358])
    ylim([56.8325 56.8405])
    pbaspect([diff(xlim)*cosd(56.8)/diff(ylim) 1 1])

end