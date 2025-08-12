addpath('../../../code/');
clear;

depname = 'F:/Alaska2022/data/iceberg_surveys/raw/20220825_singingflower/dragon/';

%tlim = parseDeploymentTimes(depname);
tlim = [datenum(2022,8,25,17,10,0) datenum(2022,8,25,18,25,0)];
adcp = parseNortekDeployment(depname,{},tlim(1),tlim(2));

%%
channel = 'burst';
beam = 1;

nt = length(adcp.(channel).time);

r_max = 2;
r_idx = adcp.(channel).range<=r_max;

%adcp.echo.amp = adcp.echo.amp - min(extrema(adcp.echo.amp));
[amp_binned,t_binned] = time_bin(adcp.(channel).time,abs(adcp.(channel).amp(:,r_idx,beam))',2/86400);
amp_binned = amp_binned'-min(extrema(amp_binned));
nt_binned = length(t_binned);

bw_binned = edge(amp_binned,'canny');

b = 2;
amp2 = b.^(10*amp_binned/max(extrema(amp_binned)));
bw2 = edge(amp2,'canny');

%% changepts
npts = 1;
stat = 'rms';

% fprintf('all data...\n')
% ipts = ones(npts,nt);
% for i = 1:nt
%     pts = findchangepts(adcp.echo.amp(i,r_idx),'statistic',stat,'maxnumchanges',npts);
%     for j = 1:length(pts)
%         ipts(j,i) = pts(j);
%     end
% end

% fprintf('binned...\n')
% ipts_binned = ones(npts,nt_binned);
% for i = 1:nt_binned
%     pts = findchangepts(amp_binned(i,:),'statistic',stat,'maxnumchanges',npts);
%     for j = 1:length(pts)
%         ipts_binned(j,i) = pts(j);
%     end
% end

fprintf('binned and amplified...\n')
ipts2 = ones(npts,nt_binned);
for i = 1:nt_binned
    pts = findchangepts(amp2(i,:),'statistic',stat,'maxnumchanges',npts);
    for j = 1:length(pts)
        ipts2(j,i) = pts(j);
    end
end

%% plot
clrs = {[0 0 0],[1 0 0],[0 1 0]};
for i = 3:3
    fprintf('plot %d...\n',i);
    figure(i); clf;
    pcolor(t_binned,adcp.(channel).range(r_idx),amp_binned')
    shading flat
    xlim(tlim)
    hold on
    datetick('x','HH:MM','keeplimits')
    colorbar
    colormap parula
    
    switch i
        case 1
            x = adcp.(channel).time;
            yi = ipts;
        case 2
            x = t_binned;
            yi = ipts_binned;
        case 3
            x = t_binned;
            yi = ipts2;
    end
    y = adcp.(channel).range(yi);
    
    % filter edge
    yi_filt = floor(meanFilter(yi(1,:),5));
    good = ~isnan(yi_filt);
    x_filt = t_binned(good);
    y_filt = adcp.(channel).range(yi_filt(good));
    
    
%     for j = 1:npts
%         plot(x,y(j,:),'.-','color',clrs{j});
%     end
    
    plot(x_filt,y_filt,'k');
end
