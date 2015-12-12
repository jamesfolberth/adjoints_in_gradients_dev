function [Wadj] = analysis_mat_adjoint_1d(lx, wname, dwname, extmode)
% analysis_mat_adjoint_1d - build wavelet analysis matrix
% 
% W = analysis_mat_adjoint_1d(lx, wname, dwname, extmode) builds the adjoint of the 
% wavelet analysis matrix for a 1d signal of length lx.
% wname is a string containing the analysis wavlet name.
% dwname is a string containing the dual wavlet name.
% extmode is a string containing the signal extension mode.
%
% Wadj is constructed via reconstruction of wavelet coefficients
% and applying the adjoint extension operator.
% 
% Supported signal extension modes are 'zpd', 'sym', and 'ppd'.
% 

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

Wadj = zeros(lx, 2*la);

for i=1:la
   ca = zeros(la,1); cd = zeros(la,1);
   ca(i) = 1;
   xe = idwt(ca, cd, dwname, lx+2*(lf-1), 'mode', 'zpd'); % take central portion
   Wadj(:,i) = extension_adjoint_1d(xe, lx, lf-1, extmode);
end
for i=la+1:2*la
   ca = zeros(la,1); cd = zeros(la,1);
   cd(i-la) = 1;
   xe = idwt(ca, cd, dwname, lx+2*(lf-1), 'mode', 'zpd'); % take central portion
   Wadj(:,i) = extension_adjoint_1d(xe, lx, lf-1, extmode);
end

end
