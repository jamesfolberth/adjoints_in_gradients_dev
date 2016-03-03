function [] = driver()

addpath('./HNO')

X = double(imread('cameraman.pgm'));
%X = double(imread('cameraman_resize.pgm'));
X = X/255; % scale to [0,1]

[P,center] = psfGauss([9,9],4);

% generate blurred image
B=imfilter(X,P,'symmetric'); 
%B=imfilter(X,P,'circular'); 
% add some Gaussian white noise
randn('seed',314);
Bobs=B + 1e-3*randn(size(B));

% show original and observed blurred image
%figure(1)
%subplot(1,2,1)
%imshow(X,[])
%title('Original')
%subplot(1,2,2)
%imshow(Bobs,[])
%title('Blur+Noise')

% FISTA parameters
%lambda = 1e-4; % used in B+T's example code
lambda = 2e-5; % used in FISTA SIAM paper

pars.MAXITER=200; % do this many iterations
pars.fig=0; % suppress the figure while running FISTA
%pars.BC='periodic';
pars.B = 1; % TODO JMF need to look this up for CDF 9/7 wavelets

% MATLAB wavelet parameters
dwtmode('sym') % symmetric (half-point) <=> reflexive
%dwtmode('per') % periodic
%dwtmode('zpd') % zero padding

% use the (bad) orthogonal DCT as the sparsifying transform
%[Xout,fun_all]=deblur_wavelet_FISTA_trans(Bobs,P,center,@dct2,@idct2,1e-4,pars);
%[Xout,fun_all]=deblur_dwt_FISTA_trans(Bobs,P,center,@wavelet_analysis,@wavelet_synthesis,lambda,pars);
[Xout,fun_all]=deblur_dwt_FISTA_trans_direct(Bobs,P,center,@wavelet_analysis,@wavelet_synthesis,@wavelet_synthesis_adjoint,lambda,pars);

% show the original and recovered images
figure(2)
subplot(1,2,1)
imshow(X,[])
title('Original')
subplot(1,2,2)
imshow(Xout,[])
title('Recovered')

fprintf(1, 'recovery l2-error (rel) = %e\n', norm(Xout-X,'fro')/norm(X,'fro'));
%fprintf(1, 'recovery nnz (%%nnz) = %d (%3.2f)\n', sum(abs(X_iter(:))>0),sum(abs(X_iter(:))>0)/numel(X_iter)*100);
%fprintf(1, 'recovery %%(big coeffs) = %3.2f\n', sum(abs(X_iter(:)) > 1e-4)/numel(X_iter)*100);


end

function [Y,S] = wavelet_analysis(X)
   % DCT
   %Y = dct2(X);
   %S = [];
   
   % 3-stage Haar, orth
   %[Y,S] = wavedec2(X, 3, 'haar');

   % 3-stage Daubechies, 4 VMs, orth
   %[Y,S] = wavedec2(X, 3, 'db4');

   % 3-stage CDF 9/7
   [Y,S] = wavedec2(X, 2, 'bior4.4');
end

function [X] = wavelet_synthesis(Y,S)
   if nargin < 2
      S = [];
   end
   %X = idct2(Y);

   %X = waverec2(Y, S, 'haar');
   
   %X = waverec2(Y, S, 'db4');

   X = waverec2(Y, S, 'bior4.4');
end

function [Y] = wavelet_synthesis_adjoint(X,S)
% Adjoint of the wavelet synthesis operator
   if nargin < 2
      S = [];
   end
   %Y = dct2(X);
   
   %Y = wavedec2(X, 3, 'haar'); % orthogonal
   
   %Y = wavedec2(X, 3, 'db4'); % orthogonal

   %TODO JMF boundary condition issues!
   %TODO where are the BC issues?!
   % Update: the operators involved here are very nearly orthonormal.
   Y = wavedec2(X, 2, 'rbio4.4');
end
