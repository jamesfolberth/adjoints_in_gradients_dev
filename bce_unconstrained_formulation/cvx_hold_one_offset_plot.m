function [] = cvx_hold_one_offset_plot(type)
% Plot recovery error (while holding h or m fixed) versus
% the our estimated lengths of the h_est, m_est

if nargin == 0
   type = 'm_known';
end

% data files
sigma = 0.005;
lambda_h = 1e-1;
lambda_m = 1e-2;
TV_h = 1e-3;
TV_m = 1e-2;
addpath('..');
addpath('../2016_01_14_rar');
%s = load('../Source Waveform - Impulsive.mat');
%S = load('../2016_01_14_rar/Impulsive Data.mat');
s = load('../Source Waveform - Gaussian, Underdetermined.mat');
S = load('../2016_01_14_rar/Underdetermined Case - Essential.mat');
%s = load('../Source Waveform - Gaussian, Overdetermined.mat');
%S = load('../2016_01_14_rar/Overdetermined Case - Essential.mat');
channel = 2;

%y = S.x(channel,:).';
h_true = S.hTrue(channel,:).';
K_true = length(h_true);

m_true = (s.s).';
N_true = length(m_true);

y = conv(h_true, m_true) + sigma*randn(N_true+K_true-1,1);

offsets = 0:5:500
conv_err = zeros(size(offsets));
m_offset_raster = zeros(numel(m_true), numel(offsets));

% if we know m
if type == 'm_known'
   for i = 1:numel(offsets)
      offset = offsets(i)
      N_est = N_true-offset;
      K_est = length(y)-N_est+1;
      m_est = m_true(1:end-offset);

      m_offset_raster(1:numel(m_est),i) = m_est(:);

      A = conv_lin_op(K_est,N_est);
      D_h = spdiags([[-ones(K_est-1,1);0] [ones(K_est-1,1);1]], [0 1], K_est, K_est);
      D_m = spdiags([[-ones(N_est-1,1);0] [ones(N_est-1,1);1]], [0 1], N_est, N_est);
      
      cvx_clear
      cvx_precision low % duality gap tends to be less than 1e-5
      cvx_begin
         variable h_est(K_est);
         %minimize 1/2*sum_square(y-A*vec(h_est*m_est')) + lambda_h*norm(h_est,1)
         %minimize 1/2*sum_square(y-A*vec(h_est*m_est')) + TV_h*norm(D_h*h_est,1)
         minimize 1/2*sum_square(y-A*vec(h_est*m_est')) + lambda_h*norm(h_est,1) + TV_h*norm(D_h*h_est,1)
      cvx_end
      
      conv_err(i) = norm(y-conv(h_est,m_est),2)/norm(y,2);
      fprintf(1, 'convolution l2-error (rel) = %e\n', conv_err(i));
      
      try % I'm doing this over SSH; might timeout
         figure(1)
         plot(offsets, conv_err)
         xlabel('m offset');
         ylabel('convolution l2-error (rel)');
      end
   end
   
   [offsets(:) conv_err(:)]

   figure(1)
   plot(offsets, conv_err)
   xlabel('m offset');
   ylabel('convolution l2-error (rel)');

%% if we know h
%TODO make this like m_known case
%elseif type == 'h_known'
%   for i = 1:numel(offsets)
%      offset = offsets(i);
%      K_est = K_true-offset
%      N_est = length(y)-K_est+1
%      h_est = h_true(1:end-offset);
%
%      A = conv_lin_op(K_est,N_est);
%      D_h = spdiags([[-ones(K_est-1,1);0] [ones(K_est-1,1);1]], [0 1], K_est, K_est);
%      D_m = spdiags([[-ones(N_est-1,1);0] [ones(N_est-1,1);1]], [0 1], N_est, N_est);
%
%      cvx_begin
%
%
%      cvx_end
%   end
end


end


function [] = cvx_hold_one_recovery()


end
