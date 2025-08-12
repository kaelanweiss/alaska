% Script to examine pitch data from adcp records during each segment.
%
% KJW
% 27 Mar 2025

clear

proc_dir = 'F:/meltstake/data/proc';
fig_dir = 'F:/meltstake/figures/adcp_attitude';

ms_tbl = loadMSInfo(26:28,'manualwindows');

[dep_nums,uidx] = unique(ms_tbl.Number);
dep_names = ms_tbl.Folder(uidx);
ndeps = length(dep_nums);

% figure
figsize = [22 10];
pad = [.07 .05 .07 .12];
shift = [0 0];

% loop through deployments
for i = 1:ndeps
    seg_nums = ms_tbl.Window(ms_tbl.Number==dep_nums(i));

    % load pck edges for boundary location
    load(fullfile(proc_dir,dep_names{i},'adv/pck_edges.mat'))
    
    % loop through segments
    for j = 1:length(seg_nums)
        % load
        load(fullfile(proc_dir,dep_names{i},sprintf('adcp%d.mat',j)))

        % fit slope to pitch data
        tfit = seconds(adcp.attitude.time - adcp.attitude.time(1));
        [b,bint,~,~,stats] = regress(adcp.attitude.pitch_ice,[ones(size(tfit)) tfit],0.05);
        
        % some useful values
        pitch_mean = mean(adcp.attitude.pitch_ice);
        pitch_std = std(adcp.attitude.pitch_ice);
        delta_pitch = b(2)*tfit(end)/2;
        pitch_fit = delta_pitch*[-1 1] + pitch_mean;
        p_rate = b(2)*3600; % deg/hr
        p_r_ci = diff(bint(2,:))/2*3600;

        % distance to ice
        row_idx = ms_tbl.Number == dep_nums(i) & ms_tbl.Window == j;
        idxt = edges.time >= ms_tbl.Start(row_idx) & edges.time <= ms_tbl.End(row_idx);
        posij = edges.pos(idxt,:);
        pos_mean = mean(posij)/10;
        
        % plot
        clear ax
        fig = figure(10*(i-1)+j); clf
        setFigureSize(fig,figsize);
        for k = 1:3
            ax(k) = axes(fig,'position',axgridpos(1,3,k,pad,shift));
            hold(ax(k),'on')
            box(ax(k),'on')
        end
        % pitch
        plot(ax(1),adcp.attitude.time,adcp.attitude.pitch_ice,'k-')
        plot(ax(1),adcp.attitude.time,hannFilter(adcp.attitude.pitch_ice,60*adcp.burst.samplerate),'r-','linewidth',1)
        plot(ax(1),adcp.attitude.time([1 end]),pitch_fit,'--','color',colors(1),'linewidth',2)
        text(ax(1),.02,.96,sprintf('pitch rate: %.2f %c %.2f deg/hr',round(p_rate,2),char(177),round(p_r_ci,2)),'units','normalized','fontsize',11)
        text(ax(1),.02,.89,sprintf('mean pitch: %.2f %c %.2f deg',round(pitch_mean,2),char(177),round(pitch_std,2)),'units','normalized','fontsize',11)
        
        % roll
        plot(ax(2),adcp.attitude.time,adcp.attitude.roll_ice,'k-')
        % heading
        plot(ax(3),adcp.attitude.time,adcp.attitude.heading,'k-')
        
        % labels
        ylabel(ax(1),'pitch [deg]')
        ylabel(ax(2),'roll [deg]')
        ylabel(ax(3),'heading [deg]')
        title(ax(2),sprintf('deployment %d segment %d',dep_nums(i),j),'fontsize',12)

        % save figure
        print(fig,fullfile(fig_dir,sprintf('adcp_prh_%02d_%02d.png',dep_nums(i),j)),'-dpng','-r350')
        
        % print pitch rate
%         fprintf('%.2f (%.2f)\n',round(p_rate,2),round(p_r_ci,2))
        
        % print mean pitch
        fprintf('%.1f (%.1f)\n',round(pitch_mean,1),round(pitch_std,1))
        
        % print ice distance
%         fprintf('%.1f (%.1f) \n',round(mean(pos_mean),1),round(std(pos_mean),1))

        a = 5;
    end

end