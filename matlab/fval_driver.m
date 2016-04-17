function [] = fval_driver()

addpath('./beck_FISTA_matlab_files/HNO')

%X = double(imread('cameraman.pgm'));
%X = double(imread('cameraman_resize.pgm'));
%X = double(imread('texmos3.pgm')); % from USC-SIPI
%X = double(imread('clock.pgm'));   % from USC-SIPI
X = double(imread('resolution.pgm')); % from USC-SIPI
X = X/255; % scale to [0,1]

[P,center] = psfGauss([9,9],4);

% generate blurred image
B=imfilter(X,P,'symmetric'); 
%B=imfilter(X,P,'circular'); 
% add some Gaussian white noise
randn('seed',314);
Bobs=B + 1e-3*randn(size(B));

% from Beck+Teboulle's deblur_...
% used for blurred image reconstruction error
[m,n]=size(Bobs);
Pbig=padPSF(P,[m,n]);
trans=@(X) 1/sqrt(prod(size(Bobs)))*fft2(X);
itrans=@(X) sqrt(prod(size(Bobs)))*ifft2(X);
% computng the eigenvalues of the blurring matrix         
Sbig=fft2(circshift(Pbig,1-center));
Btrans = trans(Bobs);


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
WSyAd_pinv = @(Y) wavelet_synthesis_adjoint_2d(Y, wname, dwname, extmode, levels, 'pinv');
WSyAd = @(Y) wavelet_synthesis_adjoint_2d(Y, wname, dwname, extmode, levels, 'adjoint');

[Xout_pinv,fun_all_pinv,X_iter_pinv]=deblur_dwt_FISTA_trans_direct(Bobs,P,center,WAn,WSy,WSyAd_pinv,lambda,pars,X);
[Xout,fun_all,X_iter]=deblur_dwt_FISTA_trans_direct(Bobs,P,center,WAn,WSy,WSyAd,lambda,pars,X); %

% run forever to get fstar
%pars_tmp = pars; pars_tmp.MAXITER = 2*pars.MAXITER;
%[Xout_star,fun_all_star,X_iter]=deblur_dwt_FISTA_trans_direct(Bobs,P,center,WAn,WSy,WSyAd,lambda,pars_tmp,X);
%fstar = min(fun_all_star(:));

% show the recovered image and some other info
imwrite(Xout_pinv, 'deblurred_images/tmp_pinv.pgm');
imwrite(Xout, 'deblurred_images/tmp.pgm');
%imwrite(Xout_star, 'deblurred_images/tmp_star.pgm');


% Plot the function values vs number of iterations
%figure(5)
%iters = 1:pars.MAXITER;
%clf();
%hold on;
%semilogy(iters, fun_all-fstar, 'LineWidth', 3);
%semilogy(iters, fun_all_pinv-fstar, '--', 'LineWidth', 3);
%set(gca, 'YScale', 'log');
%hold off;
%legend('true adjoint', 'pinv approx');
%xlabel('Iteration');
%ylabel('Objective value');

%figure(5)
%iters = 1:pars.MAXITER;
%plot(iters, fun_all_pinv, iters, fun_all, 'LineWidth', 3);
%legend('pinv approx', 'true adjoint', 'Location', 'SouthEast');
%axis([1 numel(fun_all) .6 .9]);
%xlabel('Iteration');
%ylabel('SSIM');

%save('resolution_bior4.4_sym_fval.mat', 'fun_all_pinv', 'fun_all', 'fun_all_star');
%save('resolution_bior4.4_sym_psnr.mat', 'fun_all_pinv', 'fun_all');
save('resolution_bior4.4_sym_ssim.mat', 'fun_all_pinv', 'fun_all');
%fval_plotter();

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

