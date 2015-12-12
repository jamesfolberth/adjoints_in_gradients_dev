function [] = adjoint_bc_experiment_1d(N)
   
   if nargin == 0
      N = 8;
   end

   % dwtmode
   dwtmode('zpd', 'nodisp') % we'll take care of the extension; force MATLAB to use zero padded BCs

   % extension mode
   extmode = 'sym'
   %extmode = 'ppd' % don't pick 'per'
   %extmode = 'zpd'
 
   %wn = 'haar'
   %dwn = 'haar';
   
   wn = 'db2'
   dwn = 'db2';

   %wn = 'db4'
   %dwn = 'db4';
 
   %wn = 'bior4.4'
   %dwn = 'rbio4.4';
   % CDF 9/7 has a length 9 analysis lowpass filter
   %                      7 synthesis lowpass filter
   % <=> bior 4.4 has 4 vanishing moments for the analysis highpass filter
   %                  4 vanishing moments for the synthesis highpass filter
  
   %x = 0:N-1
   %[Lo_D, Hi_D] = wfilters(wn, 'd'); % decomp filters
   %lf = length(Lo_D); % really just need length of filter
   %xe = wextend('1D', extmode, x, lf-1, 'b')
   %[a,d] = dwt(xe, wn); % wavelet transform with zero padded BC's
   
   Id = eye(N);
   [Lo_D, Hi_D] = wfilters(wn, 'd'); % decomp filters
   lf = length(Lo_D); % really just need length of filter
   xe = wextend('1D', extmode, Id(:,1), lf-1, 'b');
   [a,d] = dwt(xe, wn);
   la = numel(a); ld = numel(d);
   lenC = la + ld;
   W = zeros(lenC,N);
   Ir = eye(lenC);
   Wt = zeros(N,lenC);

   % Build wavelet analysis operator for one stage wavelet transform
   for i = 1:N
      xe = wextend('1D', extmode, Id(:,i), lf-1, 'b');
      [a,d] = dwt(xe, wn);
      W(:,i) = [a;d];
   end 
   
   % Build adjoint wavelet operator
   for i = 1:la
      ei = zeros(la,1); ei(i) = 1;
      %xe = idwt(ei, zeros(ld,1), dwn);
      xe = idwt(ei, zeros(ld,1), dwn, N+2*(lf-1));
      % this is the adjoint of the extension operator
      switch extmode
      case 'zpd',
         Wt(:,i) = xe(lf:lf+N-1);

      case 'sym',
         x = xe((lf-1)+1:(lf-1)+N);
         x(1:lf-1) = x(1:lf-1) + xe(lf-1:-1:1);
         x(end-lf+2:end) = x(end-lf+2:end) + xe(end:-1:end-lf+2);
         Wt(:,i) = x;

      otherwise
         error('bad extension mode')
      end
   end
   for i = la+1:la+ld
      ei = zeros(ld,1); ei(i-la) = 1;
      %xe = idwt(zeros(la,1), ei, dwn);
      %XXX: take the central extended length portion of the signal
      xe = idwt(zeros(la,1), ei, dwn, N+2*(lf-1));
      % this is the adjoint of the extension operator
      switch extmode
      case 'zpd',
         Wt(:,i) = xe(lf:lf+N-1);

      case 'sym',
         x = xe(lf:lf+N-1);
         x(1:lf-1) = x(1:lf-1) + xe(lf-1:-1:1);
         x(end-lf+2:end) = x(end-lf+2:end) + xe(end:-1:end-lf+2);
         Wt(:,i) = x;

      otherwise
         error('bad extension mode')
      end
   end

   %W.'
   %Wt

   %figure(); spy(W.' ~= Wt);

   %eigs(W.'*W,1)

   fprintf(1, '\\|W^T-Wt\\|_inf = %f\n', norm(W.'-Wt,'inf'))
   
   %fprintf(1, '\\|W^T-Wt\\|_inf (approx coeffs) = %f\n', norm(W(1:la,:).'-Wt(:,1:la),'inf'))
   %
   %fprintf(1, '\\|W^T-Wt\\|_inf (approx coeffs, left BC) = %f\n', norm(W(1:lf-1,:).'-Wt(:,1:lf-1),'inf'))
   %fprintf(1, '\\|W^T-Wt\\|_inf (approx coeffs, right BC) = %f\n', norm(W(la-lf+1:la,:).'-Wt(:,la-lf+1:la),'inf'))

   %fprintf(1, '\\|W^T-Wt\\|_inf (detail coeffs) = %f\n', norm(W(la+1:end,:).'-Wt(:,la+1:end),'inf'))
   %fprintf(1, '\\|W^T-Wt\\|_inf (approx coeffs, left BC) = %f\n', norm(W(1:lf-1,:).'-Wt(:,1:lf-1),'inf'))
   %fprintf(1, '\\|W^T-Wt\\|_inf (approx coeffs, right BC) = %f\n', norm(W(la:la+lf-1,:).'-Wt(:,la:la+lf-1),'inf'))



end

