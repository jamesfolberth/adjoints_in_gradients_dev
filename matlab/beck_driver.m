function [Xout] = beck_driver()

addpath('./beck_FISTA_matlab_files/HNO')

%X = double(imread('cameraman.pgm'));
%X = double(imread('cameraman_resize.pgm'));
%X = double(imread('texmos3.pgm')); % from USC-SIPI
%X = double(imread('clock.pgm'));   % from USC-SIPI
%X = double(imread('resolution.pgm')); % from USC-SIPI
X = double(imread('resolution_crop.pgm'));
X = X/255; % scale to [0,1]

[P,center] = psfGauss([9,9],4);

% generate blurred image
B=imfilter(X,P,'symmetric'); 
%B=imfilter(X,P,'circular'); 
% add some Gaussian white noise
randn('seed',314);
Bobs=B + 1e-3*randn(size(B));

% unblurred and blurred side-by-side
%blur_fig = [X(:,1:size(X,2)/2) Bobs(:,1:size(X,2)/2)];
%blur_fig(:,size(X,2)/2-3:size(X,2)/2+3) = 0;
%size(blur_fig)
%imwrite(blur_fig, 'resolution_blurred_figure.pgm');

% from Beck+Teboulle's deblur_...
% used for blurred image reconstruction error
[m,n]=size(Bobs);
Pbig=padPSF(P,[m,n]);
trans=@(X) 1/sqrt(prod(size(Bobs)))*fft2(X);
itrans=@(X) sqrt(prod(size(Bobs)))*ifft2(X);
% computng the eigenvalues of the blurring matrix         
Sbig=fft2(circshift(Pbig,1-center));
Btrans = trans(Bobs);


% show original and observed blurred image
%figure(1)
%%subplot(1,2,1)
%%imshow(X,[])
%%title('Original')
%%subplot(1,2,2)
%%imshow(Bobs,[])
%%title('Blur+Noise')
%imshow(Bobs,[])
%error('stop here')

% FISTA parameters
%lambda = 1e-4; % used in B+T's example code
lambda = 2e-5; % used in FISTA SIAM paper

pars.MAXITER=2500; % do this many iterations
pars.fig=0; % suppress the figure while running FISTA
%pars.BC='periodic';
pars.B = 1; % TODO JMF need to look this up for CDF 9/7 wavelets

% MATLAB wavelet parameters
dwtmode('zpd', 'nodisp');

%extmode = 'zpd'
extmode = 'sym'

%adjoint_mode = 'pinv';
adjoint_mode = 'adjoint';

levels = 3;

%wname  = 'db1'
%dwname = 'db1';
%wname  = 'db5'  % filter length 9
%dwname = 'db5'; 
wname  = 'bior4.4' % filter lengths 9 and 7
dwname = 'rbio4.4'; 
%wname  = 'bior2.2' % filter lengths 5 and 3
%dwname = 'rbio2.2';

L = build_wavedec_levels_2d(size(Bobs), levels, wname, extmode);

WAn = @(Y) wavelet_analysis_2d(Y, wname, extmode, levels);
WSy = @(X) wavelet_synthesis_2d(X, L, size(Bobs), wname, extmode, levels);
WSyAd = @(Y) wavelet_synthesis_adjoint_2d(Y, wname, dwname, extmode, levels, adjoint_mode);

%[Xout,fun_all]=deblur_dwt_FISTA_trans_direct(Bobs,P,center,WAn,WSy,WSyAd,lambda,pars);
[Xout,fun_all,X_iter]=deblur_dwt_FISTA_trans_direct(Bobs,P,center,WAn,WSy,WSyAd,lambda,pars,X); % X is used for SSIM

% show the original and recovered images
%figure(2)
%subplot(1,2,1)
%imshow(X,[])
%title('Original')
%subplot(1,2,2)
%imshow(Xout,[])
%title('Recovered')

fprintf(1, 'recovery l2-error (rel) = %e\n', norm(Xout-X,'fro')/norm(X,'fro'));
fprintf(1, 'blurred l2-error (rel) = %e\n', norm(Sbig.*trans(Xout) - Btrans,'fro')/norm(Btrans,'fro'));
fprintf(1, 'recovery nnz (%%nnz) = %d (%3.2f)\n', sum(abs(X_iter(:))>0),sum(abs(X_iter(:))>0)/numel(X_iter)*100);
%fprintf(1, 'recovery %%(big coeffs) = %3.2f\n', sum(abs(X_iter(:)) > 1e-4)/numel(X_iter)*100);

% show the recovered image and some other info
figure(3)
imshow(Xout,[])
imwrite(Xout, 'deblurred_images/tmp.pgm');
title(sprintf('Recovered - iter=%d, wname=''%s'', extmode=''%s''', pars.MAXITER, wname, extmode));

% Plot the decay of non-zero values in wavelet coeffs
%figure(4);
%wc = sort(abs(X_iter(:)),'descend');
%semilogy(wc);
%axis([0 1e5 1e-10 1e2]);

% Plot the function values vs number of iterations
%figure(5)
%fstar = norm(Sbig.*trans(X)-Btrans,'fro')^2+lambda*sum(sum(abs(WAn(X))));
%semilogy(1:pars.MAXITER, fun_all-fstar);
%xlabel('Iteration');
%ylabel('Objective value');

end


function [X] = wavelet_analysis_2d(Y, wname, extmode, levels)

   lY = size(Y);
   if numel(lY) == 1
      lY = [lY(1) lY(1)];
   elseif numel(lY) > 2 || numel(lY) == 0
      error('lY should have at most 2 entries')
   end
   
   assert(levels >= 1, 'Number of decomposition levels should be >= 1.');
   assert(levels <= wmaxlev(lY(1), wname), 'Number of decomposition levels too high.');
   assert(levels <= wmaxlev(lY(2), wname), 'Number of decomposition levels too high.');
   
   [Lo_D, Hi_D] = wfilters(wname, 'd'); % decomp filters
   lf = length(Lo_D);
  
   dwtmode('zpd', 'nodisp'); % wavedec2 doesn't take the 'mode' arg like (i)dwt2
   Ye = wextend('2D', extmode, Y, lf-1, 'b');
   [X,L] = wavedec2(Ye, levels, wname);

end

function [Y] = wavelet_synthesis_2d(X, L, lY, wname, extmode, levels)

   [Lo_D, Hi_D] = wfilters(wname, 'd'); % decomp filters
   lf = length(Lo_D);
  
   dwtmode('zpd', 'nodisp'); % wavedec2 doesn't take the 'mode' arg like (i)dwt2
   Ye = waverec2(X, L, wname);
   %lY = L(end,:);  % L coming from our build_wavedec_levels accounts for extension!
   Y = extension_pinv_2d(Ye, lY, lf-1, extmode);

end

function [X] = wavelet_synthesis_adjoint_2d(Y, wname, dwname, extmode, levels, adjoint_mode)

   lY = size(Y);
   if numel(lY) == 1
      lY = [lY(1) lY(1)];
   elseif numel(lY) > 2 || numel(lY) == 0
      error('lY should have at most 2 entries')
   end
   
   assert(levels >= 1, 'Number of decomposition levels should be >= 1.');
   assert(levels <= wmaxlev(lY(1), wname), 'Number of decomposition levels too high.');
   assert(levels <= wmaxlev(lY(2), wname), 'Number of decomposition levels too high.');
   
   [Lo_D, Hi_D] = wfilters(wname, 'd'); % decomp filters
   lf = length(Lo_D);
  
   dwtmode('zpd', 'nodisp'); % wavedec2 doesn't take the 'mode' arg like (i)dwt2
   if strcmp(adjoint_mode,'pinv')
      Ye = wextend('2D', extmode, Y, lf-1, 'b'); % test: this is "close" to adjoint of pinv
      [X,L] = wavedec2(Ye, levels, wname);       % test: this is "close" to adjoint of W
   elseif strcmp(adjoint_mode,'adjoint')
      Ye = extension_pinv_adjoint_2d(Y, lY, lf-1, extmode);
      [X,L] = wavedec2(Ye, levels, dwname);
   else
      error('Bad adjoint mode');
   end

end

