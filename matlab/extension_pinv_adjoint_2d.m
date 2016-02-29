function [xe] = extension_pinv_adjoint_2d(x, lX, le, extmode)
% extension_pinv_adjoint_2d - adjoint of pinv of 2d signal extension operator
% 
% xe = extension_pinv_adjoint_2d(x, lX, le, extmode) applies the adjoint of 
% pseudoinverse of the 2d  signal extension operator to the matrix xe.  lX 
% is the size of the original signal.  le is the length of the extension, 
% which is applied to both sides. extmode is a string containing the signal
% extension mode (a la wextend).

switch extmode
case 'zpd'
   xe = [zeros(lX(1),le) x zeros(lX(1),le)];
   xe = [zeros(le,lX(2)+2*le); xe; zeros(le,lX(2)+2*le)];

case 'sym'
   % this feels like a Kronecker product
   % do the first dimension extension 
   xe = zeros(lX(1)+2*le,lX(2)+2*le);
   
   dim_2_inds = le+1:le+lX(2);
   xe(1:le,                  dim_2_inds) = 0.5*x(le:-1:1,:);
   xe(le+1:2*le,             dim_2_inds) = 0.5*x(1:le,:);
   xe(2*le+1:lX(1),          dim_2_inds) = x(le+1:lX(1)-le,:);
   xe(lX(1)+1:lX(1)+le,      dim_2_inds) = 0.5*x(lX(1)-le+1:lX(1),:);
   xe(lX(1)+le+1:lX(1)+2*le, dim_2_inds) = 0.5*x(lX(1):-1:lX(1)-le+1,:);

   % clean up any overlap between extensions
   if 2*le > lX(1)
      for i=2*le:-1:lX(1)+1
         xe(2*le-i+1, dim_2_inds) = 1/3*x(i-le,:);
         xe(i,        dim_2_inds) = 1/3*x(i-le,:);
         xe(lX(1)+i,  dim_2_inds) = 1/3*x(lX(1)-i+le+1,:);
      end
   end

   % then do the second dimension extension
   xe(:, 1:le)                  = 0.5*xe(:,2*le:-1:le+1);
   xe(:, le+1:2*le)             = 0.5*xe(:,le+1:2*le);
   xe(:, 2*le+1:lX(2))          = xe(:, 2*le+1:lX(2));
   xe(:, lX(2)+le+1:lX(2)+2*le) = 0.5*xe(:, lX(2)+le:-1:lX(2)+1); % switch order to avoid overwrite
   xe(:, lX(2)+1:lX(2)+le)      = 0.5*xe(:, lX(2)+1:lX(2)+le);

   % clean up any overlap between extensions
   if 2*le > lX(2)
      % To avoid copying xe from after the first dim extension, we recomputed
      % only the necessary extended columns.
      for i=2*le:-1:lX(2)+1
         xe(:, 2*le-i+1) = 1/3*extension_pinv_adjoint_1d(x(:, i-le), lX(1), le, 'sym');
         xe(:, i) = 1/3*extension_pinv_adjoint_1d(x(:, i-le), lX(1), le, 'sym');
         xe(:, lX(2)+i) = 1/3*extension_pinv_adjoint_1d(x(:, lX(2)-i+le+1), lX(1), le, 'sym');
      end
   end
   
   %xe_check = epa_2d_sym_check(x,lX,le)
   %norm(xe-xe_check,'fro')
   
otherwise
   error('Unsupported signal extension mode: %s', extmode)
end

end

% Check symmetric extension with the following
function [Xe] = epa_2d_sym_check(X,lX,le)
   Xe = zeros(lX(1)+2*le,lX(2)+2*le);
   for j=le+1:le+lX(2)
      Xe(:,j) = extension_pinv_adjoint_1d(X(:,j-le),lX(1),le,'sym');
   end
   for i=1:lX(1)+2*le
      Xe(i,:) = extension_pinv_adjoint_1d(Xe(i,le+1:le+lX(2)),lX(2),le,'sym');
   end
end
