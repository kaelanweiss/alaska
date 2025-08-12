function instruments = parseRBRInstrumentFile(depname)
% Function to parse "instruments.txt" file in RBR folders.
% 
% instruments = parseRBRInstrumentFile(depname)
% 
% Input
%   depname: path to folder containing /RBR/ directory
% 
% Output:
%   instruments: structure with fields
%       sn: serial number
%       pos: position of instrument
%       model: RBR instrument model (solo, duet, concerto, etc)
%       t_offset: time offset for aligning clocks
%
% KJW
% 30 Aug 2022

% Open file
fid = fopen(fullfile(depname,'RBR','instruments.txt'),'r');

% Instantiate arrays
sn = [];
pos = [];
model = {};
t_offset = [];

% Read lines
while ~feof(fid)
    line = fgetl(fid);
    if line(1) ~= '#'
        C = strsplit(line,',');
        sn = [sn;str2double(C{1})];
        pos = [pos;str2double(C{2})];
        model{end+1} = C{3};
        t_offset = [t_offset;str2double(C{4})];
    end
end

% Close file
fclose(fid);

% Arrange output
instruments(length(sn)) = struct();
for i = 1:length(sn)
    instruments(i).sn = sn(i);
    instruments(i).pos = pos(i);
    instruments(i).model = model{i};
    instruments(i).t_offset = t_offset(i);
end