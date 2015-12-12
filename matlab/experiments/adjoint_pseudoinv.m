function [] = adjoint_pseudoinv()
% How close is W^*W to I?  What about WW^*?

   % dwtmode
   %dwtmode('sym')
   dwtmode('ppd')
   %dwtmode('zpd')

   N = 30;
   num_levels = 2
 
   %wn = 'haar'
   %dwn = 'haar';

   %wn = 'db4'
   %dwn = 'db4';
 
   wn = 'bior4.4'
   dwn = 'rbio4.4';
   % CDF 9/7 has a length 9 analysis lowpass filter
   %                      7 synthesis lowpass filter
   % <=> bior 4.4 has 4 vanishing moments for the analysis highpass filter
   %                  4 vanishing moments for the synthesis highpass filter
   
   Id = eye(N);
   [C,L] = wavedec(Id(:,1), num_levels, wn);
   lenC = numel(C)
   W = zeros(lenC,N);
   Ir = eye(lenC);
   Wt = zeros(N,lenC);

   % This looks like the correct procedure modulo boundary conditions,
   % which I don't yet know how to handle.  It works for zero-padding BCs.
   for i = 1:N
      [C,~] = wavedec(Id(:,i), num_levels, wn);
      W(:,i) = C;
   end 
   for i = 1:lenC
      Wt(:,i) = waverec(Ir(:,i), L, dwn);
   end
 
   %W
   %Wt.'
   %figure(); spy(W ~= Wt.');
   
   %Wt*W

   imagesc(Wt*W); colorbar() % analysis -> synthesis
   title('W^TW')
   %imagesc(W*Wt); colorbar() % synthesis -> analysis
   %title('WW^T')

   fprintf(1, '\\|W^T-Wt\\|_inf = %f\n', norm(W.'-Wt,'inf'))
   fprintf(1, '\\|I-W^TW\\|_2 = %f\n', norm(Id-Wt*W))
   fprintf(1, '\\|I-WW^T\\|_2 = %f\n', norm(Ir-W*Wt))

end
