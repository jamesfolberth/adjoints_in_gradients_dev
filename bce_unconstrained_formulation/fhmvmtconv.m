function cc = fhmvmtconv(c,r)
% cc = fhmvmtconv(c,r)
%
% Hankel to Toeplitz-ciruclant conversion
%
% Inputs
%   c   first column of the Hankel
%   r   last row of the Hankel
% Ourput
%   cc  Toeplitz-circulant 
% Note. When c(m) does not match and r(1), c(m) wins
%       over r(1).

% Kevin Browne     McMaster Univ.  November 2007

m = length(c);
n = length(r);

%
if (c(m) ~= r(1))
    disp('Col wins over row in fhmvmtconv.m')
end

% first, change Hankel to Toeplitz by reversing columns
% the first column of Toeplitz is [c(n)...c(m) r(2)...r(n)]
% and the first row of Toeplitz is [c(n) c(n-1)...c(1)]
% then the Toeplitz is expanded into a Toeplitz-circulant
% the first column of the circulant is:
%    [c(n)...c(m) r(2)...r(n) c(1)...c(n-1)], when m>=n
%    [r(n-m+1)...r(n) c(1)...c(m) r(2)...r(n-m)], otherwise

cc = zeros(m+n-1,1);             % initial first col of circulant

if (m>=n)
    for i = 1:(m-n+1)
        cc(i) = c(n+i-1);
    end
    for i = 1:(n-1)
        cc(m-n+i+1) = r(i+1);
        cc(m+i) = c(i);
    end
else % m<n
    for i = 1:m
        cc(i) = r(n-m+i);
        cc(m+i) = c(i);
    end
    for i = 1:(n-m-1)
        cc(2*m+i) = r(i+1);
    end
end

