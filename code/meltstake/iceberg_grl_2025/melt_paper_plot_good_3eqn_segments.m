% Script to plot the velocity during the 5 segments for which the
% 3-equation parameterization gets it right.
%
% KJW
% 16 Apr 2024

clear

% load data
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

% segments
s = [2 6 7 8 31];

proc_path = 'F:meltstake\data\proc';

for i = 1:length(s)
    
    load(fullfile(proc_path,ms_tbl.Folder{s(i)},sprintf('adcp%d.mat',ms_tbl.Window(s(i)))))

    adcpQuickPlot(figure(1000+i),adcp,'vel_ice',.1*[-1 1],nan,[0 1],1)

end
