%clc, close all, clear all

% Matlab script to plot a single data package from decoded ping 360 log file
% TO USE:
% 1.) download python script at
%     https://github.com/bluerobotics/ping-viewer/blob/master/examples/decode_sensor_binary_log.py
%     and install package bluerobotics-ping
% 2.) run python3 decode_sensor_binary_log.py [encoded .bin file path] > [output file path]
% 3.) specify input filepath below for decoded file.

% USER INPUT ---------------
input_filepath = 'test_ping_log.decoded';
% --------------------------

raw_data = importdata(input_filepath);


ping.num_payloads = length(find(contains(raw_data,'Payload')));
ping.data = zeros(1200,ping.num_payloads);
ping.sample_period = zeros(1,ping.num_payloads);
ping.angle = zeros(1,ping.num_payloads);
ping.timestamp = cell(1,ping.num_payloads);


for payload = 1:ping.num_payloads

    % get intensity values
    ping.data(:,payload) = get_data(raw_data, payload, 'data');
    
    % get sample_period
    ping.sample_period(payload) = get_data(raw_data, payload, 'sample_period');
    
    % get angle in degrees & remove offset s.t. arrow on ping360 is at 0 deg
    angle_grad = get_data(raw_data, payload, 'angle');
    ping.angle(payload) = 0.9*angle_grad - 200;

    % get time stamp
    ping.timestamp{payload} = get_timestamp(raw_data, payload);

    disp(payload/ping.num_payloads)
end

%%

%plot intensity vs distance
data_pt = 99;

% calculate distance in meters for each sample point
c = 1500;
sample_inc = 25*10^-9;
dist = (0.5*c*sample_inc)*ping.sample_period(data_pt)*(0:length(ping.data(:,data_pt))-1)';

set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
figure(1)
plot(dist,ping.data(:,data_pt),'-k')
title( {join(['angle = $',num2str(ping.angle(data_pt)),' ^\circ$']),join(['timestamp: $',ping.timestamp{data_pt},'$'])})
xlabel('Distance $(m)$')
ylabel('Intensity')
set(gca,'fontsize', 14) 


%%
function out = get_data(raw_data, payload, data_type)
% outputs numeric values of a specified data type for a given payload

    data_options = ["mode","gain_setting","angle","transmit_duration",...
        "sample_period","transmit_frequency","number_of_samples",...
        "data_length","data"];

    if ~any(max(strcmp(data_type,data_options)))
        error("unrecognized data type")
    end

    data_inds = find(contains(raw_data,join(['- ',data_type,': '])));

    data = raw_data(data_inds(payload));
    data = split(data{1},": ");
    out = eval(erase(data{2},"'"))';

end


function time_out = get_timestamp(raw_data, payload)
% outputs time as a string hh:mm:ss.zzz

    data_inds = find(contains(raw_data,'timestamp: '));

    data = raw_data(data_inds(payload));
    data = split(data{1},[":",".","\"]);
    times = [];
    time_ind = [3,4,7,8,11,12,15,16,17];
    for i = 1:9
        times(i) = eval(['0',data{time_ind(i)}]);
    end

    time_out = join([times(1),times(2),":",times(3),times(4),":",times(5),times(6),".",times(7),times(8),times(9)]);

end

