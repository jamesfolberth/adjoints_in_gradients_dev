%{
Simplified script to run l1 and TV
Stephen Becker, April 12 2016
%}

X   = imread('peppers.png');
X   = rgb2gray(X);
X   = double(X)/255;
N   = 40;
M   = N;
X   = X(201:200+M, 201:200+N);
imagesc(X); colormap('gray'); axis image; axis off

% Make some kind of observation matrix (e.g., a blur), for now just use
% anything simple:
m   = round( M*N/5 );
A   = randn(m,M*N);
b   = A*X(:);
%% Solve in TFOCS
addpath ~/Repos/TFOCS

% If we have a wavelet operator W and its adjoint Wt, use linop_handles
% For now, make a random matrix W
Wmatrix     = randn(N*M);

% Or ...
[C,S]       = wavedec2( X, 3, 'db2' );
Wmatrix     = zeros( length(C(:)), M*N );
E   = zeros(M,N);
for i = 1:(M*N)
    E(i)    = 1;
    [C,S]       = wavedec2( E, 3, 'db2' );
    Wmatrix(:,i)    = C(:);
    E(i)    = 0;
end
Wmatrix     = bsxfun( @times, Wmatrix, 1./sqrt(sum(Wmatrix.^2)) );

W_1         = linop_matrix(Wmatrix); % not efficient...
% celldisp( W_1(1,0) ) % size

% And make the 2D TV linear operator
W_2         = linop_TV( [M,N] );
%  W_2(1,0)  % size

% Guess some parameters
epsilon     = 1e-3; % constrains ||Ax-b||<= epsilon
alpha       = 1;
beta        = 100;
mu          = 1e0;
x0          = pinv(A)*b; % not same as A\b
% x0          = X(:); % cheating!
%%
opts        = struct('maxIts',5e4,'tol',1e-5,'debug',false,'restart',Inf);
opts.errFcn = @(fcn,dual,x) norm( A*x-b );
z0          = [];
[ x, out, opts ] = solver_sBPDN_WW( A, alpha, W_1, beta, W_2, b, epsilon, mu, x0, z0, opts );

%% print output
figure(1);
subplot(1,2,1);
imagesc(X); colormap('gray'); axis image; axis off
subplot(1,2,2);
imagesc(reshape(x,M,N)); colormap('gray'); axis image; axis off

%% compare objective values (without mu term)
alpha*norm( W_1(x,1), 1 ) + beta*norm( W_2(x,2), 1 )
alpha*norm( W_1(X(:),1), 1 ) + beta*norm( W_2(X(:),2), 1 )
