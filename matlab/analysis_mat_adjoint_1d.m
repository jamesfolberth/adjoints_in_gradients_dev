function [Wadj] = analysis_mat_adjoint_1d(lx, wname, dwname, extmode, levels)
% analysis_mat_adjoint_1d - build wavelet analysis matrix
% 
% W = analysis_mat_adjoint_1d(lx, wname, dwname, extmode) builds the adjoint of the 
% wavelet analysis matrix for a 1d signal of length lx.
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
assert(levels >= 1, 'Number of decomposition levels should be >= 1.');
assert(levels <= wmaxlev(lx, wname), 'Number of decomposition levels too high.');

[Lo_D, Hi_D] = wfilters(wname, 'd'); % decomp filters
lf = length(Lo_D);

switch extmode
case {'zpd', 'sym', 'ppd'}
   % length of approx (and detail) coeff vector for _extended_ signal, which
   % is of length lf-1 + lx + lf-1 = lx + 2(lf-1).  See dwt help page.
   la = floor((lx+3*(lf-1))/2);
otherwise
   error('Unsupported signal extension mode: %s', extmode)
end

L = build_wavedec_levels_1d(lx, levels, wname, extmode);
num_coeffs = sum(L(1:end-1)); % number of coefficients in the full decomposition
Wadj = zeros(lx, num_coeffs);

dwtmode('zpd', 'nodisp'); % waverec doesn't take the 'mode' arg like (i)dwt
for i=1:num_coeffs
   C = zeros(num_coeffs,1); C(i) = 1;
   xe = waverec(C, L, dwname);
   Wadj(:,i) = extension_adjoint_1d(xe, lx, lf-1, extmode);
end

end
