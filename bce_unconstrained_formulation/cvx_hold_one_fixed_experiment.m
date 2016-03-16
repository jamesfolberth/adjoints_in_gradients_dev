close all;
clear all;

% be sure to 
% cd ~/MATLAB/cvx
% cvx_startup

rng(0);

% random data
%K_true = 100;
%N_true = 300;
%sigma = 0.005;
%lambda_h = 1e-2;
%lambda_m = 1e-2;
%h_true = randn(K_true,1); ph = randperm(K_true); h_true(ph(1:90)) = 0;
%m_true = randn(N_true,1); pm = randperm(N_true); m_true(pm(1:270)) = 0;
%y = conv(h_true, m_true) + sigma*randn(N_true+K_true-1,1);

% data files
sigma = 0.005
lambda_h = 1e-1;
lambda_m = 1e-2;
TV_h = 1e-3;
TV_m = 1e-2;
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

y = conv(h_true, m_true) + sigma*randn(N_true+K_true-1,1);

% if we know h
%offset = 100; % >= 0
%K_est = K_true-offset
%N_est = length(y)-K_est+1
%h_est = h_true(1:end-offset);

% if ew know m
offset = 225; % >= 0
N_est = N_true-offset
K_est = length(y)-N_est+1
m_est = m_true(1:end-offset);

A = conv_lin_op(K_est,N_est);
D_h = spdiags([[-ones(K_est-1,1);0] [ones(K_est-1,1);1]], [0 1], K_est, K_est);
D_m = spdiags([[-ones(N_est-1,1);0] [ones(N_est-1,1);1]], [0 1], N_est, N_est);

cvx_begin
   %variable m_est(N_est);

   %%minimize 1/2*sum_square(y-A*vec(h_est*m_est')) + lambda_m*norm(m_est,1)
   %%minimize 1/2*sum_square(y-A*vec(h_est*m_est')) + TV_m*norm(D_m*m_est,1)
   %minimize 1/2*sum_square(y-A*vec(h_est*m_est')) + lambda_m*norm(m_est,1) + TV_m*norm(D_m*m_est,1)

   variable h_est(K_est);

   %minimize 1/2*sum_square(y-A*vec(h_est*m_est')) + lambda_h*norm(h_est,1)
   %minimize 1/2*sum_square(y-A*vec(h_est*m_est')) + TV_h*norm(D_h*h_est,1)
   minimize 1/2*sum_square(y-A*vec(h_est*m_est')) + lambda_h*norm(h_est,1) + TV_h*norm(D_h*h_est,1)

cvx_end

fprintf(1, 'convolution l2-error (rel) = %e\n', norm(y-conv(h_est,m_est),2)/norm(y,2));

figure(1)
subplot(2,1,1)
plot(y)
title('y');
subplot(2,1,2)
plot(conv(h_est,m_est));
title('conv(h\_est,m\_est)');

xmin = 1; xmax = max(numel(h_true), numel(h_est));
ymin = min(min(h_true), min(h_est)); ymax = max(max(h_true), max(h_true));
figure(2)
subplot(2,1,1)
plot(h_true)
axis([xmin xmax ymin ymax]);
title('h\_true');
subplot(2,1,2);
plot(h_est);
axis([xmin xmax ymin ymax]);
title('h\_est');

xmin = 1; xmax = max(numel(m_true), numel(m_est));
ymin = min(min(m_true), min(m_est)); ymax = max(max(m_true), max(m_true));
figure(3)
subplot(2,1,1)
plot(m_true)
axis([xmin xmax ymin ymax]);
title('m\_true');
subplot(2,1,2);
plot(m_est);
axis([xmin xmax ymin ymax]);
title('m\_est');

