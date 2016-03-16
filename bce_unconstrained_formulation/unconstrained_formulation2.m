%close all;
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

nc = 2; % number of channels

%y = S.x(channel,:).';
h_true = S.hTrue(1:nc,:).';
K_true = size(h_true,1)

m_true = (s.s).';
N_true = length(m_true)

n_zpd = 150;
h_true = [zeros(n_zpd,nc); h_true];
K_true = K_true + n_zpd;

% noise with std 0.005 is what Brendan used in his simulations
y = zeros(K_true+N_true-1, nc);
for i=1:nc
   y(:,i) = conv(h_true(:,i), m_true) + 0.005*randn(K_true+N_true-1,1);
end

%K = S.l
%N = length(y)-K+1
offset = -100;
K = K_true+offset
N = N_true-offset

% l1 terms
lambda_h = 1e-1;
lambda_m = 1e-2;
lambda = [lambda_h*ones(K*nc,1); lambda_m*ones(N,1)];
%lambda_h_vec = [zeros(100,1); 1e-2*ones(K-100,1)];
%lambda = [lambda_h_vec; lambda_m*ones(N,1)];

% TV (approx via huber) terms
lambda_h_TV = 1e-2; % too big seems to cause issues
huber_d = 0.1;

%K = 100;
%N = 300;
%sigma = 0.005;
%lambda_h = 1e-2;
%lambda_m = 1e-2;
%h_true = randn(K,1); ph = randperm(K); h_true(ph(1:90)) = 0;
%m_true = randn(N,1); pm = randperm(N); m_true(pm(1:270)) = 0;
%y = conv(h_true, m_true) + sigma*randn(N+K-1,1);

objective = @(x) unconstrained_objective2(x, nc, K, N, y, lambda_h_TV, huber_d);

%tic
%objective([randn(K,1); randn(N,1)]);
%toc

%h0 = h_true; m0 = m_true; x0 = [h0(:); m0];
h0 = randn(K,nc); m0 = randn(N,1); x0 = [h0(:); m0];
%h0 = repmat([ones(100,1); zeros(K-100,1)], [1 nc]); m0 = .5*randn(N,1); x0 = [h0(:); m0];
y0 = zeros(K+N-1,nc);
for i=1:nc
   y0(:,i) = conv(h0(:,i), m0);
end
fprintf(1, 'inital convolution l2-error (rel) = %e\n', norm(y-y0,'fro')/norm(y,'fro'));

opts = {};
opts.maxIter = 5000;
opts.optTol = 1e-2;
opts.verbose = 1;

%x_est = x0; %XXX just for testing!
%x_est = L1General2_PSSgb(objective, x0, lambda, opts);
x_est = L1General2_PSSsp(objective, x0, lambda, opts);
%x_est = L1General2_BBSG(objective, x0, lambda, opts);
h_est = reshape(x_est(1:nc*K), [K nc]);
m_est = x_est(nc*K+1:nc*K+N);

y_est = zeros(K+N-1,nc);
for i=1:nc
   y_est(:,i) = conv(h_est(:,i), m_est);
end

fprintf(1, 'inital convolution l2-error (rel) = %e\n', norm(y-y0,'fro')/norm(y,'fro'));
fprintf(1, 'convolution l2-error (rel) = %e\n', norm(y-y_est,'fro')/norm(y,'fro'));

plot_channel = 1;
figure(1);
subplot(2,1,1);
plot(h_true(:,plot_channel), 'LineWidth', 1.5);
xlim([0 max(size(h_true,1), size(h_est,1))])
title('h\_true');
subplot(2,1,2);
plot(h_est(:,plot_channel), 'LineWidth', 1.5);
xlim([0 max(size(h_true,1), size(h_est,1))])
title('h\_est');
%subplot(3,1,3);
%%plot(h_true-h_est_cheat);
%title('h\_true-h\_est\_cheat');

figure(2);
subplot(2,1,1);
plot(m_true, 'LineWidth', 1.5);
xlim([0 max(size(m_true,1), size(m_est,1))])
title('s\_true');
subplot(2,1,2);
plot(m_est, 'LineWidth', 1.5);
xlim([0 max(size(m_true,1), size(m_est,1))])
title('s\_est');
%subplot(3,1,3);
%%plot(m_true-m_est_cheat);
%title('m\_true-m\_est\_cheat');

figure(3);
subplot(2,1,1);
plot(y(:,plot_channel), 'LineWidth', 1.5);
title('x\_true');
subplot(2,1,2);
plot(y_est(:,plot_channel), 'LineWidth', 1.5);
title('x\_est = conv(h\_est, s\_est)');

