function y = fhmvmultiply(cc,x)
% y = fhmvmultiply(cc,x)
%
% Fast general Toeplitz-circulant matrix-vector multiplication
% used in fast Hankel matrix-vector multiplication
%
% Inputs
%   cc  Toeplitz-circulant converted from a Hankel
%   x   vector to be multiplied by the Hankel
% Ourput
%   y   product of the Hankel matrix and x
% Note. When c(m) does not match and r(1), c(m) wins
%       over r(1).

% Kevin Browne     McMaster Univ.  November 2007

n = length(x);
m = length(cc) - n + 1;

%
% reverse x
xrev = x(n:-1:1);

% expand xrev
xx = [xrev; zeros(m-1,1)];
%
% use fft to multiply circulant matrix-vector
yy = ifft(fft(cc).*fft(xx));

% extract the product
y = yy(1:m);
