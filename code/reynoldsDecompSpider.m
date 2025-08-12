% Script to perform a Reynolds Decomposition on Spider ADCP observations.
%
% KJW
% 21 Nov 2022

%clear

% set berg name
bergs = {'20220824_teaparty',...
         '20220824_singingflower',...
         '20220825_singingflower'};
berg = bergs{2};

% load
mat_path = 'F:/Alaska2022/data/iceberg_surveys/mat';
depname_mat = fullfile(mat_path,berg,'spider');
load(fullfile(depname_mat,'adcp.mat'))

% qa
qa_corr = adcp.burst.cor>=50;
adcp.burst.vel(~qa_corr) = nan;

%% plot deployment velocity
beam_lbls = {'right','up','left','down'};
iR = 40;
figure(1); clf;
for i = 1:4
    subplot(4,1,i)
    axlist(i) = gca;
    pcolor(adcp.burst.time,adcp.burst.range(1:iR),adcp.burst.vel(:,1:iR,i)')
    shading flat
    datetick('x','mmmdd HH:MM','keeplimits')
    cmocean('balance')
    clim(0.03*[-1 1])
    ylabel('adcp range (m)','fontsize',10)
    title(beam_lbls{i})
    hold on
end
% colorbar
pos1 = get(axlist(1),'position');
pos4 = get(axlist(4),'position');
cbar = colorbar;
cbar.Position = [pos4(1)+pos4(3)+0.005 pos4(2) 0.01 pos1(2)+pos1(4)-pos4(2)];
cbar.Label.String = ['\bf\leftarrow \rminward' repmat(' ',[1 10]) 'm/s' repmat(' ',[1 10]) 'outward \bf\rightarrow'];
cbar.Label.FontSize = 10;

% berg label
ax1pos = get(axlist(1),'position');
berg_lbl = text(axlist(1),min(xlim),max(ylim)+0.15*diff(ylim),strrep(sprintf('spider - %s',berg),'_','\_'));
set(axlist(1),'position',ax1pos);

linkaxes(axlist,'x')

%windows = [];
%% get windows (repeat as many times as necessary)
% [tget,~] = ginput(2);
% windows = [windows; tget'];

%% decomp
nD = size(windows,1);

clear chunks
chunks(nD) = struct();

for i = 1:nD
    t1 = windows(i,1);
    t2 = windows(i,2);
    idx = adcp.burst.time>=t1 & adcp.burst.time<=t2;
    nt = length(find(idx));
    chunks(i).t1 = t1;
    chunks(i).t2 = t2;
    chunks(i).tidx = idx;
    chunks(i).ice_distance = beam_distance(i,:);
    chunks(i).range = adcp.burst.range(1:iR);

    % plot time span on first figure
    for j = 1:4
        plot(axlist(j),[t1 t2],max(get(axlist(j),'ylim'))*[1 1],'k.-')
        text(axlist(j),(t1+t2)/2,max(get(axlist(j),'ylim'))-.1,num2str(i))
    end

    % mean beam velocity
    chunks(i).vel_mean = mean(adcp.burst.vel(idx,1:iR,:),1,'omitnan');

    % perturbation beam velocity
    chunks(i).vel_prime = adcp.burst.vel(idx,1:iR,:) - repmat(chunks(i).vel_mean,[nt 1 1]);
    chunks(i).vel_prime_rms = sqrt(mean(chunks(i).vel_prime.^2,1,'omitnan'));
end


%% plot
save_figs = false;
for i = 1:nD
    figure(1+i); clf;
    subplot(2,1,1); hold on
    set(gca,'clipping','off')
    for j = 1:4
        rj = chunks(i).range*min(chunks(i).ice_distance)/chunks(i).ice_distance(j);
        plot(rj,chunks(i).vel_mean(:,:,j),'.-','linewidth',1,'color',colors(mod(j,4)+1))
    end
    grid
    ylabel({'mean beam','velocity (m/s)'},'fontsize',12)
    title(sprintf('%s - section %d',strrep(berg,'_','\_'),i),'fontsize',12)
    %ylim([-.02 .02])

    subplot(2,1,2); hold on
    set(gca,'clipping','off')
    for j = 1:4
        rj = chunks(i).range*min(chunks(i).ice_distance)/chunks(i).ice_distance(j);
        plot(rj,chunks(i).vel_prime_rms(:,:,j),'.-','linewidth',1,'color',colors(mod(j,4)+1))
    end
    grid
    legend({'right','up','left','down'},'fontsize',10);
    xlabel('scaled range (m)','fontsize',12)
    ylabel({'rms beam velocity','perturbation (m/s)'},'fontsize',12)

    if save_figs
        print(figure(i+1),sprintf('../figures/weekly/20221122/%s_RD_%d.png',berg,i),'-dpng','-r400')
        set(axlist(1),'xlim',[chunks(i).t1 chunks(i).t2])
        berg_lbl.Position(1) = chunks(i).t1;
        for j = 1:4
            datetick(axlist(j),'x','mmmdd HH:MMZ','keeplimits')
            set(axlist(j),'clim',0.03*[-1 1])
        end
        print(figure(1),sprintf('../figures/weekly/20221122/%s_TS_%d.png',berg,i),'-dpng','-r400')
    end

end



