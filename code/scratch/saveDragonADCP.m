% script to quickly parse Dragon CTD and ADCP
% could trim them down in time
%
% KJW
% 2023/10/25

deps = {'20230921_180214','20230922_163549','20230925_193312'};
dir_raw = 'G:\Shared drives\Ice-ocean-interactions\fieldwork_docs_and_data\Leconte2309\data\raw\dragon';


tlims = {{'20230921_2030','20230921_2330'},{'20230922_1730','20230922_2330'},{'20230925_2015','20230925_2145'}};

for i = 1:length(deps)
    t1 = datenum(tlims{i}{1},'yyyymmdd_HHMM');
    t2 = datenum(tlims{i}{2},'yyyymmdd_HHMM');
    adcp = parseNortekDeployment(fullfile(dir_raw,deps{i},'adcp'),t1,t2);
    save(fullfile(dir_raw,deps{i},'adcp','adcp'),'adcp')
    figure(i); clf
    plot(adcp.attitude.time,adcp.attitude.pressure)
    datetick('x','mmmdd HH:MM','keeplimits')
end