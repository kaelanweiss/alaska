% Script to reparse 881a meltstake deployments
%
% KJW
% 3 Feb 2026
clear

raw_dir = 'F:/meltstake/data/raw';

dep_tbl = loadMSInfo(20:25);
ndeps = size(dep_tbl,1);

for i = 1:ndeps
    depname = dep_tbl.Folder{i};
    t1 = dep_tbl.Start(i)-minutes(30);
    t2 = dep_tbl.End(i)+minutes(30);
    fprintf('%s\n',depname)
    try
        sonar = parse881aDeployment(fullfile(raw_dir,depname),t1,t2);
    catch
        fprintf('\tfailed to find/parse data\n')
        continue
    end

    % save
    fprintf('\tsaving...\n')
    save(fullfile(raw_dir,depname,'sonar881a','sonar881a.mat'),'sonar')

end
