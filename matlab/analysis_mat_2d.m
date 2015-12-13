function [W] = analysis_mat_2d(lX, wname, extmode, levels)
% analysis_mat_2d - build wavelet analysis matrix
% 
% W = analysis_mat_2d(lX, wname, extmode) builds  the wavelet analysis matrix
% for a rank-2 array of size lX.  wname is a string containing the wavlet name.
% extmode is a string containing the signal extension mode.
%
% W is constructed by computing the wavelet coefficients for the ith canonical
% basis vector e_i.  Each column of W contains coefficients organized like the 
% output of wavedec.
% 
% Supported signal extension modes are 'zpd', 'sym', and 'ppd'.  The signal is
% extended on all sides.
%
% TODO: this could maybe be extended to wavelet packets, maybe?
% 

if nargin < 4
   levels = 1;
end

if numel(lX) == 1
   lX = [lX(1) lX(1)];
elseif numel(lX) > 2 || numel(lX) == 0
   error('lX should have at most 2 entries')
end

assert(levels >= 1, 'Number of decomposition levels should be >= 1.');
assert(levels <= wmaxlev(lX(1), wname), 'Number of decomposition levels too high.');
assert(levels <= wmaxlev(lX(2), wname), 'Number of decomposition levels too high.');

[Lo_D, Hi_D] = wfilters(wname, 'd'); % decomp filters
lf = length(Lo_D);

switch extmode
case {'zpd', 'sym', 'ppd'}
   % length of approx (and detail) coeff vector for _extended_ signal, which
   % is of length lf-1 + lX(i) + lf-1 = lX(i) + 2(lf-1).  See dwt2 help page.
   lA = floor((lX+3*(lf-1))/2);
otherwise
   error('Unsupported signal extension mode: %s', extmode)
end

% compute total number of coefficients (used for pre-allocation)
nrows = 3*prod(lA);
lA_l = lA;
for l=2:levels
   lA_l = floor((lA_l+(lf-1))/2);
   nrows = nrows + 3*prod(lA_l);
end
nrows = nrows + prod(lA_l); % final approx coeffs

W = zeros(nrows, prod(lX));

dwtmode('zpd', 'nodisp'); % wavedec2 doesn't take the 'mode' arg like (i)dwt2
for i=1:prod(lX);
   x = zeros(lX); x(i) = 1; % put a 1 at the ith entry of the matrix
   xe = wextend('2D', extmode, x, lf-1, 'b');
   
   [C,L] = wavedec2(xe, levels, wname);
   W(:,i) = C(:); % unroll the whole thing
end

end
