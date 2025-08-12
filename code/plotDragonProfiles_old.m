% Script to plot profiles from a Dragon deployment. This is also a
% prototype for how to process the data streams.
%
% KJW
% 8 Sep 2022
clear;

raw_path = 'F:/Alaska2022/data/iceberg_surveys/raw';
berg = '20220825_singingflower';

depname = fullfile(raw_path,berg,'dragon');

% load profile file
times = load(fullfile(depname,'profile_times.mat'));
n = length(times.t1);
t_start = times.t1(1)-60/86400;
t_end = times.t2(end)+60/86400;

% RBR
fprintf('RBR...\n');
rbr = parseRBRDeployment(depname,t_start,t_end);
rake_idx = ~strcmp('concerto',extractfield(rbr,'model'));
ctd_idx = find(~rake_idx,1);
pos = extractfield(rbr(rake_idx),'pos');

% ROV
fprintf('ROV...\n');
rov = parseROVTelemetry(depname,t_start,t_end);

% ADCP
fprintf('ADCP...\n');
adcp = parseNortekDeployment(depname,{},t_start,t_end);
if ~isempty(fieldnames(adcp))
    use_adcp = true;
else
    use_adcp = false;
end


%% build rake profiles (this should be transitioned into a processing script)
% time on vertical axis
clear profs
profs(n) = struct();
inst_idx = find(rake_idx);
t_pad = .015/86400;
for i = 1:n
    idx = rbr(inst_idx(1)).time>=times.t1(i)-t_pad & rbr(inst_idx(1)).time<=times.t2(i)+t_pad;
    profs(i).time = rbr(inst_idx(1)).time(idx);
    % ctd
    profs(i).ctd = struct();
    profs(i).ctd.depth = interp1(rbr(ctd_idx).time,rbr(ctd_idx).values(:,strcmp(rbr(ctd_idx).channels,'Depth')),profs(i).time);
    profs(i).ctd.temperature = interp1(rbr(ctd_idx).time,rbr(ctd_idx).values(:,strcmp(rbr(ctd_idx).channels,'Temperature')),profs(i).time);
    profs(i).ctd.salinity = interp1(rbr(ctd_idx).time,rbr(ctd_idx).values(:,strcmp(rbr(ctd_idx).channels,'Salinity')),profs(i).time);
    % rov
    profs(i).rov = struct();
    profs(i).rov.heading = interp1angle(rov.time,rov.heading*pi/180,profs(i).time)*180/pi;
    profs(i).rov.pitch = interp1angle(rov.time,rov.pitch*pi/180,profs(i).time)*180/pi;
    profs(i).rov.roll = interp1angle(rov.time,rov.roll*pi/180,profs(i).time)*180/pi;
    
    % clean up heading to minimize wrapping
    hdg_mode = mode(round(profs(i).rov.heading));
    profs(i).rov.heading = profs(i).rov.heading - hdg_mode;
    if max(abs(diff(profs(i).rov.heading)))>350
        cuts = min(extrema(profs(i).rov.heading)):max(extrema(profs(i).rov.heading));
        performance = NaN*(cuts);
        for j = 1:length(cuts)
            new_hdg = profs(i).rov.heading;
            if cuts(j) < 0
                new_hdg(new_hdg<cuts(j)) = new_hdg(new_hdg<cuts(j))+360;
            else
                new_hdg(new_hdg>=cuts(j)) = new_hdg(new_hdg>=cuts(j))-360;
            end
            performance(j) = max(diff(new_hdg));
        end
        [~,cut_idx] = min(performance);
        cut_best = cuts(cut_idx);
        if cut_best < 0 
            profs(i).rov.heading(profs(i).rov.heading<cut_best) = profs(i).rov.heading(profs(i).rov.heading < cut_best)+360;
        else
            profs(i).rov.heading(profs(i).rov.heading>=cut_best) = profs(i).rov.heading(profs(i).rov.heading>=cut_best)-360;
        end
    end
    
    % remove mode from pitch and roll
    profs(i).rov.pitch = profs(i).rov.pitch - mode(round(profs(i).rov.pitch,1));
    profs(i).rov.roll = profs(i).rov.roll - mode(round(profs(i).rov.roll,1));
    
    
    % rake
    m = length(profs(i).time);
    profs(i).T = NaN(m,length(pos));
    profs(i).T(:,1) = rbr(inst_idx(1)).values(idx,strcmp(rbr(inst_idx(1)).channels,'Temperature'));
    for j = 2:length(inst_idx)
        idx = rbr(inst_idx(j)).time>=times.t1(i)-t_pad & rbr(inst_idx(j)).time<=times.t2(i)+t_pad;
        profs(i).T(:,j) = rbr(inst_idx(j)).values(idx,strcmp(rbr(inst_idx(j)).channels,'Temperature'));
    end
end

%% plot
% berg title
berg_splt = strsplit(berg,'_');
berg_ttl = sprintf('%s %s (local)',berg_splt{2},datestr(datenum(berg_splt{1},'yyyymmdd'),'dd-mmm-yyyy'));
% time axis
lw = 1.2;
fs = 14;
for i = [2 3 4 11 12]
    t_sec = 86400*(profs(i).time-profs(i).time(1));
    %t_sec = profs(i).ctd.depth;
    figure(i); clf;
    ax = [];
    % temp
    subplot(1,7,1:3);
    ax(end+1) = gca;
    pcolor([pos 2*pos(end)-pos(end-1)] ,t_sec,[profs(i).T profs(i).T(:,end)]);
    shading flat;
    axis ij;
    ax1pos = get(gca,'position');
    cbar = colorbar;
    cbar.Position = [sum(ax1pos([1 3]))+.0025 ax1pos(2) 0.01 ax1pos(4)];
    ylabel('profile time (s)','fontsize',fs);
    xlabel('rake position (cm)','fontsize',fs)
    title([berg_ttl ' |  profile ' num2str(i)],'fontsize',fs+2);
    %caxis([4.8 9]);
    cmocean('thermal');    
    
    % depth
    subplot(1,7,4); hold on;
    ax(end+1) = gca;
    plot(profs(i).ctd.depth,t_sec,'k.-','linewidth',lw);
    idx = diff(profs(i).ctd.depth)<=0;
    plot(profs(i).ctd.depth(idx),t_sec(idx),'ro');
    axis ij;
    grid;
    box on;
    set(gca,'YTickLabel',{});
    %set(gca,'XTick',[0:5:max(profs(i).ctd.depth)+5]);
    xlabel('depth (m)','fontsize',fs);
    
    % ctd
    S_offset = round(min(profs(i).ctd.salinity)-min(profs(i).ctd.temperature));
    subplot(1,7,5);
    ax(end+1) = gca;
    hold on;
    plot(profs(i).ctd.temperature,t_sec,'linewidth',lw);
    plot(profs(i).ctd.salinity-S_offset,t_sec,'linewidth',lw);
    axis ij;
    grid;
    box on;
    set(gca,'YTickLabel',{});
    %set(gca,'XTick',[0:5:max(profs(i).depth)+5]);
    xlabel('CTD (psu,C)','fontsize',fs);
    legend({'T',['S-' num2str(S_offset)]},'fontsize',fs-2,'location','best');
    xlim(extrema(union(profs(i).ctd.temperature,profs(i).ctd.salinity-S_offset))+[-0.25 0.25]);
    
    % heading
    subplot(1,7,6);
    ax(end+1) = gca;
    plot(profs(i).rov.heading,t_sec,'k-','linewidth',lw);
    axis ij;
    grid;
    set(gca,'YTickLabel',{});
    xlabel('heading (deg)','fontsize',fs);
    xlim(extrema(profs(i).rov.heading));
    
    % roll/pitch
    rp_thresh = 5;
    subplot(1,7,7);
    ax(end+1) = gca;
    hold on;
    set(gca,'YColor','k');
    plot(profs(i).rov.pitch,t_sec,'k-','linewidth',lw);
    plot(profs(i).rov.pitch(abs(profs(i).rov.pitch)>rp_thresh),t_sec(abs(profs(i).rov.pitch)>rp_thresh),'r.','linewidth',lw);
    plot(profs(i).rov.roll,t_sec,'k--','linewidth',lw);
    plot(profs(i).rov.roll(abs(profs(i).rov.roll)>rp_thresh),t_sec(abs(profs(i).rov.roll)>rp_thresh),'r.','linewidth',lw);
    axis ij;
    grid;
    box on;
    set(gca,'YAxisLocation','right')
    xlabel('pitch/roll (deg)','fontsize',fs);
    legend({'pitch','roll'},'fontsize',fs-2,'location','southwest');
    
    linkaxes(ax,'y');
end

%%
% dur = 5;
% for i = 1:n
%     figure(i);
%     fprintf('%d) ',i);
%     for j = 1:dur
%         fprintf('%d...',dur+1-j)
%         pause(1)
%     end
%     fprintf('0\n')
% end