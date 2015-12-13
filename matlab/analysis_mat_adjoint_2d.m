function [Wadj] = analysis_mat_adjoint_2d(lX, wname, dwname, extmode, levels)
% analysis_mat_adjoint_2d - build wavelet analysis matrix
% 
% W = analysis_mat_adjoint_2d(lx, wname, dwname, extmode) builds the adjoint of the 
% wavelet analysis matrix for a 2d signal of length lx.
% wname is a string containing the analysis wavlet name.
% dwname is a string containing the dual wavlet name.
% extmode is a string containing the signal extension mode.
%
% Wadj is constructed via dual reconstruction of wavelet coefficients
% and applying the adjoint extension operator.
% 
% Supported signal extension modes are 'zpd', 'sym', and 'ppd'.
% 

if nargin < 5
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

L = build_wavedec_levels_2d(lX, levels, wname, extmode);
% number of coefficients in the full decomposition
num_coeffs = prod(L(1,:)) + 3*sum(prod(L(2:levels+1,:),2));
Wadj = zeros(prod(lX), num_coeffs);

dwtmode('zpd', 'nodisp'); % waverec2 doesn't take the 'mode' arg like (i)dwt2
for i=1:num_coeffs
   C = zeros(num_coeffs,1); C(i) = 1;
   xe = waverec2(C, L, dwname);
   x = extension_adjoint_2d(xe, lX, lf-1, extmode);
   Wadj(:,i) = x(:);
end

end
