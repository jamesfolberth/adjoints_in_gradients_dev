close all;
clear all;

%rng(0);

% Mark Schmidt's 2012 minFunc
%addpath('./minFunc2012/minFunc');
%addpath('./minFunc2012/autoDif');

% Mark Schmidt's L1General
cd 'L1General'
addpath(genpath(pwd));
cd ..

% data files
addpath('..');
addpath('../2016_01_14_rar');
s = load('../Source Waveform - Impulsive.mat');
S = load('../2016_01_14_rar/Impulsive Data.mat');
%s = load('../Source Waveform - Gaussian, Underdetermined.mat');
%S = load('../2016_01_14_rar/Underdetermined Case - Essential.mat');
%s = load('../Source Waveform - Gaussian, Overdetermined.mat');
%S = load('../2016_01_14_rar/Overdetermined Case - Essential.mat');
channel = 2;

%y = S.x(channel,:).';
h_true = S.hTrue(channel,:).';
K_true = length(h_true)

m_true = (s.s).';
N_true = length(m_true)

% noise with std 0.005 is what Brendan used in his simulations
y = conv(h_true, m_true) + 0.005*randn(N_true+K_true-1,1);

%K = S.l
%N = length(y)-K+1
offset = 2;
K = K_true+offset;
N = N_true-offset;


lambda_h = 1e-1; % we know plausible channels
lambda_h_vec = [zeros(100,1); 1e-2*ones(K-100,1)];
lambda_m = 1e-1;
lambda = [lambda_h*ones(K,1); lambda_m*ones(N,1)];
%lambda = [lambda_h_vec; lambda_m*ones(N,1)];

%K = 100;
%N = 300;
%sigma = 0.005;
%lambda_h = 1e-2;
%lambda_m = 1e-2;
%h_true = randn(K,1); ph = randperm(K); h_true(ph(1:90)) = 0;
%m_true = randn(N,1); pm = randperm(N); m_true(pm(1:270)) = 0;
%y = conv(h_true, m_true) + sigma*randn(N+K-1,1);

objective = @(x) unconstrained_objective(x, K, N, y);

%tic
%objective([randn(K,1); randn(N,1)]);
%toc

%h0 = randn(K,1); m0 = randn(N,1); x0 = [h0; m0];
h0 = [ones(100,1); zeros(K-100,1)]; m0 = .5*randn(N,1); x0 = [h0; m0]; % seems to help
%h0 = [ones(100,1); zeros(K-100,1)]; m0 = 1e-1*ones(N,1); x0 = [h0; m0];
%x0 = [h_true; m_true] + 0.5*randn(K+N,1);
%x0 = [h_true; m_true] + 0.0*randn(K+N,1);
% Use 1st singular vectors to get good initial guess for l2 term
%A = conv_lin_op(K,N);
%tmp = reshape(A\y, [K N]); %TODO eventually want to do CG or something to find this, not building sparse mat
%[U,S,V] = svd(tmp, 'econ');
%h0 = sqrt(S(1,1))*U(:,1);
%m0 = sqrt(S(1,1))*V(:,1);
%x0 = [sqrt(S(1,1))*U(:,1); sqrt(S(1,1))*V(:,1)];
fprintf(1, 'inital convolution l2-error (rel) = %e\n', norm(y-conv(h0,m0),2)/norm(y));

opts = {};
opts.maxIter = 5000;
opts.verbose = 1;
%x_est = x0; %XXX just for testing!
x_est = L1General2_PSSsp(objective, x0, lambda, opts);
h_est = x_est(1:K);
m_est = x_est(K+1:K+N);

% CHEATING! Let's try to get the right constant
%I = find(h_est ~= 0 & h_true ~= 0); % find nonzero shared by h_est and h_true
%if numel(I) > 0
%   h_est_cheat = h_est./h_est(I(1))*h_true(I(1));
%else
%   h_est_cheat = h_est;
%   fprintf(1, 'too sparse to cheat!\n');
%end
%
%I = find(m_est ~= 0 & m_true ~= 0); % find nonzero shared by m_est and m_true
%if numel(I) > 0
%   m_est_cheat = m_est./m_est(I(1))*m_true(I(1));
%else
%   m_est_cheat = m_est;
%   fprintf(1, 'too sparse to cheat!\n');
%end

fprintf(1, 'inital convolution l2-error (rel) = %e\n', norm(y-conv(h0,m0),2)/norm(y));
fprintf(1, 'convolution l2-error (rel) = %e\n', norm(y-conv(h_est,m_est),2)/norm(y,2));
%fprintf(1, 'h_est_cheat l1-error (rel) = %e\n', norm(h_true-h_est_cheat,1)/norm(h_true,1));
%fprintf(1, 'm_est_cheat l1-error (rel) = %e\n', norm(m_true-m_est_cheat,1)/norm(m_true,1));

figure(1);
subplot(2,1,1);
plot(h_true);
title('h\_true');
subplot(2,1,2);
plot(h_est);
title('h\_est');
%subplot(3,1,3);
%%plot(h_true-h_est_cheat);
%title('h\_true-h\_est\_cheat');

figure(2);
subplot(2,1,1);
plot(m_true);
title('s\_true');
subplot(2,1,2);
plot(m_est);
title('s\_est');
%subplot(3,1,3);
%%plot(m_true-m_est_cheat);
%title('m\_true-m\_est\_cheat');

figure(3);
subplot(2,1,1);
plot(y);
title('x\_true');
subplot(2,1,2);
plot(conv(h_est,m_est));
title('x\_est = conv(h\_est, s\_est)');

