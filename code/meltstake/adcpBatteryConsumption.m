clear

%% setup
data_path = 'F:/AK_202307/adcp';

deployments = {[1],...
               [1 2 3 4],...
               [1 2 3 4 5]};

vbat(2) = struct('sn',[],'t',[],'V',[]);
vbat(1).sn = 102620;
vbat(2).sn = 103041;
sn_list = [vbat.sn];

for i = 1:length(deployments)
    ndeps = length(deployments{i});

    for j = 1:ndeps
        adcp = parseNortekDeployment(fullfile(data_path,sprintf('ms%02d/%d',i,deployments{i}(j))),'none');
        snj = adcp.cfg.sn;

        sn_idx = sn_list==snj;

        vbat(sn_idx).t = [vbat(sn_idx).t; adcp.attitude.time; nan];
        vbat(sn_idx).V = [vbat(sn_idx).V; adcp.attitude.batteryvoltage; nan];
    end
end

for i = 1:length(vbat)
    [vbat(i).t,idx_sort] = unique(vbat(i).t);
    vbat(i).V = vbat(i).V(idx_sort);
end

%%
figure(100); clf; hold on
for i = 1:2
    plot(vbat(i).t,hannFilter(vbat(i).V,25))
end