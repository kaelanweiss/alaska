function out = hannFilter(in,k,varargin)
% Function: Hanning window filter of width k for 1D array in. If 'omitnan'
% is chosen, the each output value is normalized according to the number of
% available data points within the Hanning window of each output value.
%
% out = hannFilter(in,k)
% out = hannFilter(in,k,nanflag)
%
% Input
%   in: vector to be filtered
%   k: width of window
%   nanflag: (optional) 'omitnan' to ignore nans in mean calculation
%
% Output
%   out: filtered vector

% find nans
idx_nan = isnan(in);

% set nans to zero in original vector
y0 = in;
y0(idx_nan) = 0;

% smooth using convolution of normalized Hanning window
k = k + 2; % hanning window has zeros on both sides, don't want to count those
w = hann(k)/sum(hann(k));
y0s = conv(y0,w,'same');

% calculate normalization by convolving logical ~nan vector
y1 = 1*~idx_nan;
y1s = conv(y1,w,'same');

% if there are too many nans, keep the nan value
nan_limit = 0.3;
y1s(y1s<nan_limit) = nan;

% normalize convolved input
out = y0s./y1s;