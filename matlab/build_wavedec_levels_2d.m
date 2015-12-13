function [L] = build_wavedec_levels_2d(lX, levels, wname, extmode)
% build_wavedec_levels_2d - build the levels structure L used by wavedec2/waverec2
% 
% TODO
% This accounts for our extension of the signal
% Would be easy to extend to m x n x 3 images
%

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

L = zeros(levels+2, 2);
L(levels+2,:) = lX(:) + 2*(lf-1);
L(levels+1,:) = lA(:);

lA_l = lA;
for l=2:levels
   lA_l = floor((lA_l+(lf-1))/2);
   L(levels+2-l,:) = lA_l(:);
end
L(1,:) = lA_l(:);

end
