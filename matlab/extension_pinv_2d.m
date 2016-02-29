function [x] = extension_pinv_2d(xe, lX, le, extmode)
% extension_pinv_2d - pseudoinverse of 2d signal extension operator
% 
% x = extension_pinv_2d(xe, lX, le, extmode) applies the pseudoinverse of the 2d
% signal extension operator to the matrix xe.  lX is the size of the original
% signal.  le is the length of the extension, which is applied to both sides.
% extmode is a string containing the signal extension mode (a la wextend).
%

switch extmode
case 'zpd'
   x = xe(le+1:le+lX(1), le+1:le+lX(2));

case 'sym'
   % this feels like a Kronecker product
   % do the first dimension update (fold and add, like in 1d)
   x = xe(le+1:le+lX(1), :);
   x(1:le, :) = 0.5*(x(1:le, :) + xe(le:-1:1, :));
   x(end-le+1:end, :) = 0.5*(x(end-le+1:end, :) + xe(end:-1:end-le+1, :));
   
   % then do the second dimension update (fold and add)
   xe = x; % alias
   x = x(:, le+1:le+lX(2));
   x(:, 1:le) = 0.5*(x(:, 1:le) + xe(:, le:-1:1));
   x(:, end-le+1:end) = 0.5*(x(:, end-le+1:end) + xe(:, end:-1:end-le+1));
   
otherwise
   error('Unsupported signal extension mode: %s', extmode)
end

end
