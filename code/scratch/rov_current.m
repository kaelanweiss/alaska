raw_dir = 'F:/AK_202309/data/unsorted_rov/video';

ass_files = dir(fullfile(raw_dir,'*.ass'));
nfiles = length(ass_files);

% fname = '2023-09-23_19.14.43.ass';
% fname = '2023-09-23_18.46.02.ass';
% fname = '2023-09-22_00.06.59.ass';
% fname = '2023-09-23_19.14.43.ass';

delim = '{\\pos(1615,1075)}';
current = nan(10000,1);
pos = 1;

for i = 1:nfiles
    fname = ass_files(i).name;
    fprintf('(%d/%d) %s\n',i,nfiles,fname)

    fid = fopen(fullfile(raw_dir,fname),'r');
    while ~feof(fid)
        line = fgetl(fid);
        C1 = strsplit(line,delim);
        if length(C1)==2
            C2 = strsplit(C1{2},' A');
            current(pos) = str2double(C2{1});
            pos = pos + 1;
        end
    end
end

%%
figure(1); clf
plot(current,'.-')

figure(2); clf
histogram(current(current>8),0:50)

figure(3); clf
histogram(current,0:50)