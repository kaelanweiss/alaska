clear;

bergs = {'20220822_oyster','20220825_singingflower'};
berg_num = 2;
berg = bergs{berg_num};

switch berg_num
    case 1
        adcp_dp = 0.45;
    case 2
        adcp_dp = 0.35;
end

load(fullfile('F:/Alaska2022/data/iceberg_surveys/proc',berg,'dragon','profiles.mat'));

%% plot some adcp stuff
pn = 9;
beam_dirs = {'right/down','right/up','left/up','left/down'};
flds = {'vel','cor','amp'};
units = {'m/s','%','dB'};
cmaps = {'balance','matter','amp'};
clims = {0.12*[-1 1],[20 100],[53 87]};
bt_color = {0*[1 1 1],1*[1 1 1],1*[1 1 1]};

ttl_bits = strsplit(berg,'_');
berg_ttl = sprintf('%s %s',ttl_bits{2},datestr(datenum(ttl_bits{1},'yyyymmdd'),'dd-mmm'));
fs = 14;
for i = 1:3
    figure(i); clf;
    % adcp beams
    for j = 1:profs(pn).adcp.burst.nbeams
        subplot(1,profs(pn).adcp.burst.nbeams,j); hold on;
        if strcmp(flds{i},'xxxxx') % remove BT
            tmp = profs(pn).adcp.burst.(flds{i})(:,:,j)-...
                repmat(interp1(profs(pn).adcp.bt.time,profs(pn).adcp.bt.vel(:,j),profs(pn).adcp.burst.time,'nearest','extrap'),[1 profs(pn).adcp.burst.ncells]);
        else
            tmp = profs(pn).adcp.burst.(flds{i})(:,:,j);
        end
        pcolor(profs(pn).adcp.burst.range,profs(pn).adcp.burst.pressure-adcp_dp,tmp);
        plot(profs(pn).adcp.bt.distance(:,j),profs(pn).adcp.bt.pressure-adcp_dp,'.','color',bt_color{i});        
        
        shading flat
        cmocean(cmaps{i})
        %colorbar
        xlabel('range (m)','fontsize',fs);
        if j==1
            ylabel('depth (m)','fontsize',fs);
        end
        axis ij
        title(beam_dirs{j},'fontsize',fs);
        
        xlim(extrema(profs(pn).adcp.burst.range)+[0 2]);
        xlim([0 2])
        ylim(extrema(profs(pn).adcp.burst.pressure-adcp_dp));
        %set(gca,'xdir','reverse');
        caxis(clims{i});
        box on;
    end
    % add colorbar to last panel
    ax_end = gca;
    cbar = colorbar;    
    
    % echo sounder
%     subplot(1,profs(pn).adcp.burst.nbeams+1,j+1);
%     pcolor(profs(pn).adcp.echo.range,profs(pn).adcp.echo.pressure-adcp_dp,profs(pn).adcp.echo.amp)
%     shading flat
%     axis ij
%     colorbar
%     xlim([0 4])
%     caxis([50 110])
%     title('echo','fontsize',fs);
%     xlabel('range (m)','fontsize',fs)
%     %ylabel('depth (m)','fontsize',fs)
%     set(gca,'xdir','reverse');
    
    %suptitle(sprintf('%s profile %d %s',berg_ttl,pn,flds{i}));
    %suptitle(sprintf('%s: profile %d %s',flds{i},pn,berg_ttl));
    
    cbar.Position = cbarpos(ax_end,.0008,.0125);
    cbar.Label.String = units{i};
    cbar.Label.FontSize = fs;
%     cbar.Label.Rotation = 0;
%     cbar.Label.Position = [.5 1.05 0];
end
figure(1);
%%
att_flds = {'heading','pitch','roll','ahrsgyro','batteryvoltage'};
units = {'deg','deg','deg','deg/s','V'};
lw = 1.2;
figure(4); clf;
for i = 1:length(att_flds)
    subplot(1,length(att_flds),i);
    plot(profs(pn).adcp.attitude.(att_flds{i}),profs(pn).adcp.attitude.pressure-adcp_dp,'linewidth',lw);
    ylim(extrema(profs(pn).adcp.attitude.pressure-adcp_dp));
    title(att_flds{i},'fontsize',fs);
    ylabel('depth','fontsize',fs);
    xlabel(units{i},'fontsize',fs)
    grid;
    axis ij
    if strcmp(att_flds{i},'ahrsgyro')
        legend({'x','y','z'},'fontsize',fs-2,'location','best')
    end
end
suptitle(sprintf('motion: %s profile %d',berg_ttl,pn));

%%
figure(5); clf;
hold on;
for i = 1:profs(pn).adcp.bt.nbeams
    plot(profs(pn).adcp.bt.vel(:,i),profs(pn).adcp.bt.pressure-adcp_dp);
end
axis ij
xlim(0.1*[-1 1])
ylim(extrema(profs(pn).adcp.bt.pressure-adcp_dp))
legend(beam_dirs,'fontsize',fs)
xlabel('bt velocity (m/s)','fontsize',fs)
ylabel('depth (m)','fontsize',fs)
grid;
box on;
title(sprintf('%s profile %d',berg_ttl,pn),'fontsize',fs);

figure(3);
figure(2);
figure(1);
%%
figure(6); clf;
pcolor(profs(pn).adcp.echo.range,profs(pn).adcp.echo.pressure,profs(pn).adcp.echo.amp)
shading flat
axis ij
cbar = colorbar;
xlim([0 4])
caxis([50 110])
title('echo','fontsize',fs);
xlabel('range (m)','fontsize',fs)
ylabel('depth (m)','fontsize',fs)

