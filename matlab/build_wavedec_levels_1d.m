function [L] = build_wavedec_levels_1d(lx, levels, wname, extmode)
% build_wavedec_levels_1d - build the levels structure L used by wavedec/waverec
% 
% TODO
% This accounts for our extension of the signal
%

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

L = zeros(levels+2,1);
L(levels+2) = lx+2*(lf-1);
L(levels+1) = la;

la_l = la;
for l=2:levels
   la_l = floor((la_l+(lf-1))/2);
   L(levels+2-l) = la_l;
end
L(1) = la_l;

end
