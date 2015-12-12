function [W] = analysis_mat_1d(lx, wname, extmode, levels)
% analysis_mat_1d - build wavelet analysis matrix
% 
% W = analysis_mat_1d(lx, wname, extmode) builds  the wavelet analysis matrix
% for a 1d signal of length lx.  wname is a string containing the wavlet name.
% extmode is a string containing the signal extension mode.
%
% W is constructed by computing the wavelet coefficients for the ith canonical
% basis vector e_i.  Each column of W contains [ca1; cd1; ca2; cd2; ...].
% 
% Supported signal extension modes are 'zpd', 'sym', and 'ppd'.  The signal is
% extended on both sides.
% 

if nargin < 4
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

if lx < lf
   error('Signal length is shorter than filter length.');
end

% compute total number of coefficients (used for pre-allocation)
la_l = la;
nrows = 2*la;
for l=2:levels
   la_l = floor((la_l+(lf-1))/2);
   nrows = nrows + 2*la_l;
end

W = zeros(nrows, lx);

for i=1:lx;
   x = zeros(lx,1); x(i) = 1; % ith canonical basis vector
   xe = wextend('1D', extmode, x, lf-1, 'b');
   [ca, cd] = dwt(xe, wname, 'mode', 'zpd');
   W(1:2*la,i) = [ca; cd];
   
   % handle the remaining levels
   % we don't do any extension here, since we're working off of approx coeffs
   row_pointer = 2*la;
   la_l = la;
   for l=2:levels
      la_l = floor((la_l+(lf-1))/2);
      [ca_l, cd_l] = dwt(ca, wname, 'mode', 'zpd');
      W(row_pointer+1:row_pointer+2*la_l,i) = [ca_l; cd_l];
      row_pointer = row_pointer + 2*la_l;
   end
end

end
