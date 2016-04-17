function [X_out,fun_all,X_iter]=deblur_dwt_FISTA_trans_direct(Bobs,P,center,WAn,WSy,WSyAd,lambda,pars,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements FISTA for solving the linear inverse problem with 
% an orthogonal l1 wavelet regularizer and either reflexive or periodic boundary
% conditions
%
% Based on the paper
% Amir Beck and Marc Teboulle, "A Fast Iterative Shrinkage-Threshold Algorithm
% for Linear Inverse Problems",  to appear in SIAM Journal on Imaging
% Sciences
% -----------------------------------------------------------------------
% Copyright (2008): Amir Beck and Marc Teboulle
% 
% FISTA is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
% 
% INPUT
%
% Bobs............................. The observed image which is blurred and noisy
% P .................................... PSF of the blurring operator
% center ......................  A vector of length 2 containing the center
%                                           of the PSF
%XXX JMF hacked the args
% WAn - wavelet analysis
% WSy - wavelet synthesis
% WSyAd - wavelet synthesis adjoint
% lambda ...................... Regularization parameter
% pars.................................Parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.fig ............................... 1 if the image is shown at each
%                                                      iteration, 0 otherwise (Default=1)
% pars.BC .................................. boundary conditions.
%                                                      'reflexive'  (default)  or 'periodic
% pars.B - upper frame inequality bound (B=1 for orthogonal transformation)
% 
% X .................................. Original, unblurred image.  Used only for SSIM
%
% OUTPUT
% 
% X_out ......................... Solution of the problem
%                                          min{||A(X)-Bobs||^2+lambda \|Wx\|_1
% fun_all .................... Array containing all function values
%                                          obtained in the FISTA method


% Assigning parameters according to pars and/or default values
flag=exist('pars');
if (flag&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
if(flag&isfield(pars,'fig'))
    fig=pars.fig;
else
    fig=1;
end
if (flag&isfield(pars,'BC'))
    BC=pars.BC;
else
    BC='reflexive';
end
if (flag&isfield(pars,'B'))
    B=pars.B;
else
    B=1;
end
if (flag&isfield(pars,'tv'))
    tv=pars.tv;
else
    tv='iso';
end

% If there are two output arguments, initalize the function values vector.
if (nargout>=2)
    fun_all=[];
end

ssim_K = [0.01 0.03];
ssim_window = fspecial('gaussian', 11, 1.5);
ssim_L = 1;

[m,n]=size(Bobs);
Pbig=padPSF(P,[m,n]);

switch BC
    case 'reflexive'
        trans=@(X)dct2(X);
        itrans=@(X)idct2(X);
        % computng the eigenvalues of the blurring matrix         
        e1=zeros(m,n);
        e1(1,1)=1;
        Sbig=dct2(dctshift(Pbig,center))./dct2(e1);
    case 'periodic'
        trans=@(X) 1/sqrt(m*n)*fft2(X);
        itrans=@(X) sqrt(m*n)*ifft2(X);
        % computng the eigenvalues of the blurring matrix         
        Sbig=fft2(circshift(Pbig,1-center));
    otherwise
        error('Invalid boundary conditions should be reflexive or periodic');
end
% computing the two dimensional transform of Bobs
Btrans=trans(Bobs);

%The Lipschitz constant of the gradient of ||A(X)-Bobs||^2
%TODO JMF: I winged this bound; don't know if it's right.
L=B^2*2*max(max(abs(Sbig).^2));

% initialization
X_iter = WAn(Bobs);

Y=X_iter;
t_new=1;
fprintf('************\n');
fprintf('**FISTA**\n');
fprintf('************\n');
fprintf('#iter  fun-val         relative-dif\n==============================\n');
for i=1:MAXITER
    % Store the old value of the iterate and the t-constant
   X_old=X_iter;
   t_old=t_new;
   
    % Gradient step
    WY = WSy(Y);
    D=Sbig.*trans(WY)-Btrans;
    Y=Y-2/L*WSyAd(real(itrans(conj(Sbig).*D)));
   
    % soft thresholding
    D = abs(Y)-lambda/L;
    X_iter = sign(Y).*((D>0).*D);
     
    %updating t and Y
    t_new=(1+sqrt(1+4*t_old^2))/2;
    Y=X_iter+(t_old-1)/t_new*(X_iter-X_old);
    
    % Compute the l1 norm of the wavelet transform and the function value and store it in
    % the function values vector fun_all if exists.
    t=sum(sum(abs(X_iter)));
    WX_iter = WSy(X_iter);

    %fun_val=norm(Sbig.*trans(WX_iter)-Btrans,'fro')^2+lambda*t;
    [fun_val,~] = ssim_index(X, reshape(WSy(X_iter), size(X)), ssim_K, ssim_window, ssim_L);
    %fun_val = 10*log10(prod(size(X))*1^2/norm(X-WSy(X_iter),'fro')^2);
    if (nargout>=2)
        fun_all=[fun_all;fun_val];
    end
    % printing the information of the current iteration
    fprintf('%3d    %15.5f                %15.5f \n',i,fun_val,norm(X_iter-X_old,'fro')/norm(X_old,'fro'));
    
    if (fig)
        figure(314)
        imshow(X_iter,[])
    end
end

X_out=WSy(X_iter);
