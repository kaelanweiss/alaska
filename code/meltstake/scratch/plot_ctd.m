%folder = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\LeConte2305\data\concerto';
folder = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\Leconte2307\data\rbr';

d = dir(fullfile(folder,'*95_*.rsk'));

files = {d.name}';
good = 1:length(files);

figure(101); clf; hold on
axis ij
for i = 1:length(files)
    fprintf('Loading %s...\n',files{i})
    try
        rsk = RSKopen(fullfile(folder,files{i}));
        rsk = RSKreaddata(rsk);
    catch
        fprintf('Failed to open file\n\n')
        good(i) = nan;
        continue
    end
    
    plot(dn2dt(rsk.data.tstamp),rsk.data.values(:,3)+i)
    
    fprintf('\n')
end

good(isnan(good)) = [];

legend(files(good),'location','eastoutside')