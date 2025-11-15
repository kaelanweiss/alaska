% Script to process ADV deployments in the glider tank 12 Jan 2024

angle = [0 10 20];
pck(length(angle)) = struct;

data_path = '../../../../instrumentation/adv/data/pck_tank_test';

% load data
for i = 1:3
    load(fullfile(data_path,sprintf('pck%02d',angle(i)),'adv.mat'),'adv');
    pck(i).pck_time = adv.pck_time;
    pck(i).pck_dist = adv.pck_dist;
    pck(i).pck_amp = adv.pck_amp;
    pck(i).hdr_time = adv.hdr_time;
    pck(i).dp = adv.dp;
    pck(i).dprobe_ave = adv.dprobe_ave;
    pck(i).dsvol_ave = adv.dsvol_ave;
    pck(i).angle = angle(i);
end

%% distance from transducer to wall
m = 2.5/10; % known melt rate (mm/s)
nt = 90;
dt = 10;
time_wall = dt*(0:nt-1)';
d0 = [80 80 85]; % need another offset due to transducer rotating toward wall

for i = 1:3
    pck(i).time_wall = time_wall;
    pck(i).dist_wall = (d0(i) + m*time_wall)/cosd(angle(i));
end

%% calculate distance based on pck profiles
x0 = [75 40 40];
L = 250;
L0 = 100;
for i = 1:3
    out = findADVIceEdge(pck(i),0,x0(i),L,L0,3);
    pck(i).calc_time = out.time;
    pck(i).calc_pos = out.pos;
    pck(i).calc_desc = out.desc;
    pck(i).calc_hann_width = out.hann_width;
end

%% plot my known wall distance against ADV's reported distances
% dsvol_ave
figure(1001); clf
for i = 1:3
    ti = seconds(pck(i).hdr_time - pck(i).hdr_time(1));
    subplot(1,3,i)
    hold on
    plot(ti,pck(i).dsvol_ave(:,1),'k.-')
    plot(ti,pck(i).dsvol_ave(:,2),'ko--','markersize',3)
    plot(pck(i).time_wall,pck(i).dist_wall-100,'r.-')
    title(sprintf('dist to sample vol (%ddeg)',angle(i)))
    xlabel('time (s)')
    ylabel('dist (mm)')
    box on
end

% dprobe_ave
figure(1002); clf
for i = 1:3
    ti = seconds(pck(i).hdr_time - pck(i).hdr_time(1));
    subplot(1,3,i)
    hold on
    plot(ti,pck(i).dprobe_ave(:,1),'k.')
    plot(ti,pck(i).dprobe_ave(:,2),'k.')%,'markersize',3)
    plot(pck(i).time_wall,pck(i).dist_wall+50,'r-')
    %plot(pck(i).time_wall,pck(i).dist_wall*cosd(angle(i))+50,'r--')
    tcalc = seconds(pck(i).calc_time - pck(i).calc_time(1));
    for j = 1:3
        plot(tcalc,pck(i).calc_pos(:,j),'-','color',colors(j),'linewidth',1)
    end
    title(sprintf('dist to probe (%ddeg)',angle(i)))
    xlabel('time (s)')
    ylabel('dist (mm)')
    box on
end

% dp
lgd_lbl = {'b1 (down)','b2 (near)','b3 (far)'};
figure(1003); clf
for i = 1:3
    ti = seconds(pck(i).hdr_time - pck(i).hdr_time(1));
    subplot(1,3,i)
    hold on
    clear h
    for j = 1:3
        h(j) = plot(ti,interp1(1:length(pck(i).pck_dist(1,:)),pck(i).pck_dist(1,:),pck(i).dp(:,j,1)),'k.-','color',colors(j+4));
        %plot(ti,interp1(1:length(pck(i).pck_dist(1,:)),pck(i).pck_dist(1,:),pck(i).dp(:,j,2)),'ko--','markersize',3,'color',colors(j+4))
        plot(ti,count2dist_nortek(pck(i).dp(:,j,2),1430),'ko--','markersize',3,'color',colors(j+4))
    end
    legend(h,lgd_lbl,'location','northeast')
    plot(pck(i).time_wall,pck(i).dist_wall,'r.-')
    title(sprintf('dist to probe (%ddeg)',angle(i)))
    xlabel('time (s)')
    ylabel('dist (mm)')
    box on
    %ylim([-2 250])
end

%% mother of all plots
figure(999); clf

bm_lbl = {'b1 (down)','b2 (near)','b3 (far)'};

clear ax
for i = 1:3
    for j = 1:3
        ax(3*(i-1)+j) = axes(figure(999),'position',axgridpos(3,3,3*(j-1)+i,[.05,.02,.08,.08]));
        box on
        hold on
        
        % pcolor
        t0 = pck(i).pck_time(1);
        t = seconds(pck(i).pck_time-t0);
        tmax = max(t);
        pcolor(t,pck(i).pck_dist(2,:),pck(i).pck_amp(:,:,j)')
        shading flat
        
        % calculated
        t = seconds(pck(i).calc_time-t0);
        plot(t,pck(i).calc_pos(:,j),'r-','linewidth',1)
        
        % dprobe_ave
        t = seconds(pck(i).hdr_time-t0);
        plot(t,pck(i).dprobe_ave(:,1),'w.')
        
        % actual
        t = pck(i).time_wall+30;
        plot(t,pck(i).dist_wall,'y--','linewidth',1)
        
        xlim([0 tmax])
        ylim([0 349])
        
        set(gca,'xtick',0:120:900)
        
        if j==1
            title(sprintf('angle: %d deg | m=%.3f mm/s',angle(i),m/cosd(angle(i))))
        end
        if i==1
            ylabel({bm_lbl{j},'dist (mm)'})
        else
            set(gca,'yticklabel',{})
        end
        if j==3
            xlabel('time (s)')
        else
            set(gca,'xticklabel',{})
        end
    end
end
        
        
        