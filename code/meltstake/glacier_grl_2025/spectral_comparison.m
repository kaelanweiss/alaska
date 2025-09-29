% Script to compare velocity spectra between meltstakes and available
% mooring/other velocity measurements.
%
% KJW
% 26 Sep 2025
clear

addpath('..')

% collect velocity data
% meltstake
raw_dir = 'F:/meltstake/data/raw';
ms_tbl = loadMSInfo(26:28,'manualwindows');
[dep_nums,uidx] = unique(ms_tbl.Number);
dep_names = ms_tbl.Folder(uidx);
ndeps = length(dep_nums);

% 2024
% workhorse (N1)
wh24_file = 'Y:\proc\mooring\workhorse\adcp\SN14158_54.nc';
wh24 = load_nc(wh24_file);

% 2025
% workhorse
wh25_file = 'Z:\proc\mooring\workhorse\SN15645.nc';
wh25 = load_nc(wh25_file);

function data = load_nc(file)
    finfo = ncinfo(file);
    varnames = {finfo.Variables.Name}';
    data = struct();
    for i = 1:length(varnames)
        data.(varnames{i}) = ncread(file,varnames{i});
    end


end