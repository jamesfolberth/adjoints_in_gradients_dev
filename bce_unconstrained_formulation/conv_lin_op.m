function [A] = conv_lin_op(K,N)
% Build a sparse matrix that implements the linear operator in Ahmed et al. 2013
% for signals of length K and N, respectively, and with trivial subspace models
% B = I, C = I.
%
% h is of size [K 1]
% m is of size [N 1]
%
% y = \mathcal{A}(h*m') <==> y = A*vec(h*m')
%
% We'll see how fast this is... the matrix is pretty sparse.
% I don't know what the adjoint transform is in terms of FFT things, but that's
% got to be faster.
%

L = K+N-1; % length of convolved signal

m = L; n = K*N; % size of matrix
I = []; J = []; V = [];

% down first column of h*m'
for i=1:K
   for j=1:min(i,N)
      % at index (i-j+1,j) of h*m'
      %[i-j+1 j i K*(i-j+1-1)+j]
      I = [I; i];
      J = [J; K*(j-1)+i-j+1];
      V = [V; 1];
   end
   %[]
end

% across last row of h*m'
for j=2:N
   for i=K:-1:1
      % at index (i, K-i+j) of h*m'
      if K-i+j > N
         break
      end
      %[i K-i+j]
      I = [I; K+j-1];
      J = [J; K*(K-i+j-1) + i];
      V = [V; 1];
   end
   %[]
end

A = sparse(I,J,V, m,n);

end
