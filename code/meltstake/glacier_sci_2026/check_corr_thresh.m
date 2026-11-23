% Script to quickly check that the correlation threshold set at 50 isn't
% changing the U_inst results.
clear

dep_tbl = loadMSInfo(26:28);
seg_tbl = loadMSInfo(26:28,'segments');

U50(3) = struct();
U80(3) = struct();

% calculate
for i = 1:3
    % load
    seg_tbl = loadMSInfo(dep_tbl.Number(i),'segments');
    load(fullfile('F:/meltstake/data/raw',dep_tbl.Folder{i},'adcp/adcp.mat'))

    rmin = min(seg_tbl.rmax);

    % 50% correlation
    U50(i).time = adcp.burst.time;
    U50(i).Uinst = msVelocityMax(adcp,rmin);

    % 80% correlation
    adcp = msADCPTransform(adcp,80,0);
    U80(i).time = adcp.burst.time;
    U80(i).Uinst = msVelocityMax(adcp,rmin);

    % trim
    idxt = adcp.burst.time >= seg_tbl.Start(1) & adcp.burst.time <= seg_tbl.End(end);
    U50(i).Uinst(~idxt) = nan;
    U80(i).Uinst(~idxt) = nan;

end

%% segment means
for i = 1:3
    seg_tbl = loadMSInfo(dep_tbl.Number(i),'segments');
    n = size(seg_tbl,1);
    U50(i).Useg = nan(n,1);
    U80(i).Useg = nan(n,1);
    for j = 1:n
        t1 = seg_tbl.Start(j);
        t2 = seg_tbl.End(j);
        idxt = U50(i).time >= t1 & U50(i).time <= t2;

        U50(i).Useg(j) = mean(U50(i).Uinst(idxt),'omitnan');
        U80(i).Useg(j) = mean(U80(i).Uinst(idxt),'omitnan');
    end

end

%% plot
for i = 1:3
    figure(10*i+1); clf
    hold on
    plot(U50(i).time,U50(i).Uinst)
    plot(U80(i).time,U80(i).Uinst)

    figure(10*i+2); clf
    hold on
    scatter(U50(i).Uinst,U80(i).Uinst,4,'k','filled','o')
    plot([0 1],[0 1],'r-')
    scatter(U50(i).Useg,U80(i).Useg,16,'b','filled','o')
    xlim(0.8*[0 1])
    ylim(0.8*[0 1])
end

%% fit
for i = 1:3
    fprintf('---Deploy %d---\n',i)
    x = U50(i).Useg;
    y = U80(i).Useg;
    idx_nan = isnan(x) | isnan(y);
    x(idx_nan) = [];
    y(idx_nan) = [];
    [f,g] = fit(x,y,'poly1')
end