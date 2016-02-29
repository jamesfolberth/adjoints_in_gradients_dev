function [xe] = extension_pinv_adjoint_1d(x, lx, le, extmode)
% extension_pinv_adjoint_1d - adjoint of pinv of 1d signal extension operator
% 
% xe = extension_pinv_adjoint_1d(x, lx, le, extmode) applies the adjoint of the
% pseudoinverse of the 1d signal extension operator to the vector xe.  lx is 
% the length of the original signal.  le is the length of the extension, which
% is applied to both sides.  extmode is a string containing the signal extension
% mode (a la wextend).

switch extmode
case 'zpd'
   if size(x,1) == 1 % row vec
      xe = [zeros(1,le) x zeros(1,le)];
   else % col vec
      xe = [zeros(le,1); x; zeros(le,1)];
   end

case 'sym'
   if size(x,1) == 1 % row vec
      xe = zeros(1,lx+2*le);
   else
      xe = zeros(lx+2*le,1);
   end

   xe(1:le) = 0.5*x(le:-1:1);
   xe(le+1:2*le) = 0.5*x(1:le);
   xe(2*le+1:lx) = x(le+1:lx-le);
   xe(lx+1:lx+le) = 0.5*x(lx-le+1:lx);
   xe(lx+le+1:lx+2*le) = 0.5*x(lx:-1:lx-le+1);

   % clean up any overlap between extensions
   if 2*le > lx
      for i=2*le:-1:lx+1
         %[2*le-i+1 i lx+i] % print some indexes for debugging
         %[lx+i lx-i+le+1]
         xe(2*le-i+1) = 1/3*x(i-le);
         xe(i) = 1/3*x(i-le);
         xe(lx+i) = 1/3*x(lx-i+le+1);
      end
   end

otherwise
   error('Unsupported signal extension mode: %s', extmode)
end

end

% Can test symmetric case (which is not ON) with the following
%  lx=10; le=6; x = (1:lx)'; E = [circshift(fliplr(eye(le,lx)), le,2); eye(lx); fliplr(eye(le,lx))];
%  x, xep=extension_pinv_adjoint_1d(x, lx, le, 'sym'), pinv(E)'*x, norm(xep - pinv(E)'*x)

