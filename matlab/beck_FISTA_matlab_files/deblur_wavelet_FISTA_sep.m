function [X_out,fun_all]=deblur_wavelet_FISTA_sep(Bobs,P,center,W,WT,lambda,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements FISTA for solving the linear inverse problem with 
% an orthogonal l1 wavelet regularizer and a seperable PSF
%
% Based on the paper
% Amir Beck and Marc Teboulle, "A Fast Iterative Shrinkage-Threshold Algorithm
% for Linear Inverse Problems",  to appear in SIAM Journal on Imaging
% Sciences (2008)
% ------------------------------------------------------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% 
% INPUT
%
% Bobs............................. The observed image which is blurred and noisy
% P .................................... PSF of the blurring operator
% center ......................  A vector of length 2 containing the center
%                                           of the PSF
% W ....................................  A function handle. For an image
%                                           X, W(X)  is  an orthogonal
%                                           transform of the image X.
% WT .................................  A function handle. For an image
%                                           X, WT(X) is the inverse (with respect to the operator W) 
%                                           orthogonal transform of the image X 
% lambda ...................... Regularization parameter
% pars.................................Parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.fig ............................... 1 if the image is shown at each
%                                                      iteration, 0 otherwise (Default=1)
% pars.BC .................................. boundary conditions.
%                                                      'reflexive' (default) , 'periodic or 'zero'
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
if (flag&isfield(pars,'tv'))
    tv=pars.tv;
else
    tv='iso';
end

% If there are two output arguments, initalize the function values vector.
if (nargout==2)
    fun_all=[];
end

[m,n]=size(Bobs);
Pbig=padPSF(P,[m,n]);

[m,n]=size(Bobs);
Pbig=padPSF(P,[m,n]);
% computing the terms of the kronecker product
[Ar,Ac]=kronDecomp(Pbig,center,BC);
% computing the singular values of the Kronecker product of Ar and Ac
sr=svd(Ar);
sc=svd(Ac);
Sbig=sc*sr';

%The Lipschitz constant of the gradient of ||A(X)-Bobs||^2
L=2*max(max(abs(Sbig).^2));

% initialization
X_iter=Bobs;
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
     Y=Y-2/L*Ac'*(Ac*Y*Ar'-Bobs)*Ar;
    
    % Wavelet transform
    WY=W(Y);
    % Soft thresholding 
    D=abs(WY)-lambda/(L);
    WY=sign(WY).*((D>0).*D);
    % The new iterate inverse wavelet transform of WY
    X_iter=WT(WY);
    
    
    %updating t and Y
    t_new=(1+sqrt(1+4*t_old^2))/2;
    Y=X_iter+(t_old-1)/t_new*(X_iter-X_old);
    
    % Compute the l1 norm of the wavelet transform and the function value and store it in
    % the function values vector fun_all if exists.
    t=sum(sum(abs(W(X_iter))));
    fun_val=norm(Ac*X_iter*Ar'-Bobs,'fro')^2+lambda*t;
    if (nargout==2)
        fun_all=[fun_all;fun_val];
    end
    % printing the information of the current iteration
    fprintf('%3d    %15.5f                %15.5f \n',i,fun_val,norm(X_iter-X_old,'fro')/norm(X_old,'fro'));
    
    if (fig)
        figure(314)
        imshow(X_iter,[])
    end
end

X_out=X_iter;
