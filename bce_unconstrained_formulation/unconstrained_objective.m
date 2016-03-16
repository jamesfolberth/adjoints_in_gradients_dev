function [f, df] = unconstrained_objective(x, lh, lm, y)
% x = [h; m] 
% min_{h,m} 1/2 ||y - A(h*m')||_2^2 + lambda_h ||h||_1 + lambda_m ||m||_1
%
% Here we compute the 1/2||...||_2^2 term and its gradient.
%
% I think this is restricted to h,m real.
%
% A(h*m') == conv(h,m) computes the convolution of h and m from their outer product.
%
% A*(y) takes y and makes a Hankel matrix.  TODO Not sure where that comes from 
% algebraically, but it's empirically true.  We can do a Hankel mat-vec in O(n*log(n)) 
% time with the FFT (Hankel -> Toeplitz -> embed in circulant, which is diagonalized
% by the DFT).

h = x(1:lh); m = x(lh+1:lh+lm);

Ahm = conv(h,m);
%ha = [h; zeros(lm-1,1)];
%ma = [m; zeros(lh-1,1)];
%Ahm = real(ifft(fft(ha).*fft(ma)));

% objective value
f = 1/2*norm(y-Ahm, 2)^2;

% gradients
resid = Ahm-y;

%R = reshape(At*resid, [lh lm]);
%dfh = R*m;
%dfm = R.'*h;

% compute action of adjoint
hc = resid(1:lh); hr = resid(lh:lh+lm-1); % 1st col, last row of Hankel matrix
tc = resid(lm:lh+lm-1); tr = resid(lm:-1:1);  % 1st col, 1st row of corresponding Toeplitz matrix

c = [tc; tr(end:-1:2)]; % 1st col of Toeplitz T embedded in circulant matrix

Fc = fft(c);
dfh = real(ifft(Fc.*fft([m(end:-1:1); zeros(lh-1,1)]))); dfh = dfh(1:lh);
dfm = real(ifft(conj(Fc).*fft([h; zeros(lm-1,1)]))); dfm = dfm(lm:-1:1);

df = [dfh; dfm];

end

%R = hankel(hc, hr);
%dfh = R*m;
%dfm = R.'*h;

%T = toeplitz(tc, tr); % T == R*fliplr(eye(lm)), P == fliplr(eye(lm))
%dfh = T*m(end:-1:1); % keep track of H*m == (H*P)*P'*m == T*P'*m
%dfm = T.'*h; dfm = dfm(end:-1:1);

