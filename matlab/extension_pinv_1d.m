function [x] = extension_pinv_1d(xe, lx, le, extmode)
% extension_pinv_1d - pseudoinverse of 1d signal extension operator
% 
% x = extension_pinv_1d(xe, lx, le, extmode) applies the pseudoinverse of the 1d
% signal extension operator to the vector xe.  lx is the length of the original
% signal.  le is the length of the extension, which is applied to both sides.
% extmode is a string containing the signal extension mode (a la wextend).
%

switch extmode
case 'zpd'
   x = xe(le+1:le+lx);

case 'sym'
   x = xe(le+1:le+lx);
   x(1:le) = 0.5*(x(1:le) + xe(le:-1:1));
   x(end-le+1:end) = 0.5*(x(end-le+1:end) + xe(end:-1:end-le+1));
 
otherwise
   error('Unsupported signal extension mode: %s', extmode)
end

end
