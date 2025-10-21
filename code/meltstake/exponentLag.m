function y_lag = exponentLag(y,dt,T_efold,cutoff)
% NEED TO FILL IN DOCUMENTATION

% e-folding time constant
lambda = 1/T_efold;

% number of points in convolution
N = ceil(-log(cutoff)/(dt*lambda));
if ~mod(N,2) % if even, add one so convolution is centered
    N = N+1;
end

% convolving function
f = zeros(2*N-1,1);
f(1:N) = exp((1:N)*dt*lambda);
f = flip(f/sum(f));

% convolve
y_lag = conv(y,f,'same');
y_lag(1:N) = nan;