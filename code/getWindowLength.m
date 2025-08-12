function [L,n] = getWindowLength(T,T1,T2)
% Function to find the number and length of half-overlapping windows given
% a record length (T) and lower and upper bounds (T1,T2) on window length.
% There may be multiple ways to split a record and satisfy the window
% bounds. This function returns the largest possible window. A solution is
% guaranteed if T1 <= (2/3)T2. I haven't proved this, but I think it's
% true.
%
% [L,n] = getWindowLength(T,T1,T2)
%
% Input
%   T: record length
%   T1: lower bound on window 
%   T2: upper bound on window length
%
% Output
%   L: length of window
%   n: number of windows
%
% KJW
% 30 Jan 2024

L = nan;
n = nan;

% check that T1 < T2
if T1 >= T2
    warning('T1 must be less than T2')
    return
end

% special case of T < T1
if T < T1
    L = T;
    n = 1;
    return
end

% main calculation
for nj = floor(2*T/T2-1):ceil(2*T/T1-1) % possible n values
    Lj = 2*T./(nj+1); % window length with half-overlap
    if (Lj>=T1) && (Lj<=T2) % check if window length fits
        L = Lj;
        n = nj;
        return
    end
end
