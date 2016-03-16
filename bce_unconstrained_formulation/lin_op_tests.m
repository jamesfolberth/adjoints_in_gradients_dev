clear all
rng(0);

K = 4; N = 6; L = K+N-1;
%K = randi(100); N = randi(100); L = K+N-1; % just for testing
%[K N]

A = conv_lin_op(K,N);

% Test action of A
%h = randn(K,1);
%s = randn(N,1);
%
%hs = h*s'
%X = vec(h*s');
%
%x = conv(h,s)
%A*vec(hs)
%fprintf(1, 'lin op error = %e\n', norm(x-A*X,'inf'));


% Test action of A'
x = randn(K+N-1,1)
reshape(A'*x, [K N])
x1 = x(1:K); x2 = x(K:K+N-1);
hankel(x1,x2)

