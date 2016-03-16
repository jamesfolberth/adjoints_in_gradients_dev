function [f, df] = unconstrained_objective2(x, nc, lh, lm, y, lambda_h_TV, huber_d)
% x = [h1; h2; m] 
% min_{h,m} 1/2 ||y - A(h*m')||_2^2 + lambda_h ||h||_1 + lambda_m ||m||_1
%
% Here we compute the 1/2||...||_2^2 term and its gradient.
%
% This is restricted to h_i, m real, but maybe it's possible to extend to complex signals.
%
% A(h*m') == conv(h,m) computes the convolution of h and m from their outer product.
%
% A*(y) takes y and makes a Hankel matrix.  We can do a Hankel mat-vec in O(n*log(n)) 
% time with the FFT (Hankel -> Toeplitz -> embed in circulant, which is diagonalized
% by the DFT).
%
% TODO: TV penalty (approx with huber)

h = reshape(x(1:nc*lh), [lh nc]);
m = x(nc*lh+1:nc*lh+lm);

Ahm = zeros(lh+lm-1, nc);
for i=1:nc
   Ahm(:,i) = conv(h(:,i),m);
end

% objective value
f = 1/2*norm(y-Ahm, 'fro')^2;

% gradients
resid = Ahm-y;

%R = reshape(At*resid, [lh lm]);
%dfh = R*m;
%dfm = R.'*h;

% compute action of adjoint
hc = resid(1:lh,:); hr = resid(lh:lh+lm-1,:); % 1st col, last row of Hankel matrix
tc = resid(lm:lh+lm-1,:); tr = resid(lm:-1:1,:);  % 1st col, 1st row of corresponding Toeplitz matrix

c = [tc; tr(end:-1:2,:)]; % 1st col of Toeplitz T embedded in circulant matrix

Fc = fft(c); % DFT over columns
Fma = fft([m(end:-1:1); zeros(lh-1,1)]); % DFT of extended m
Fha = fft([h; zeros(lm-1,nc)]);          % DFT of extended h

dfh = real(ifft(bsxfun(@times, Fc, Fma))); dfh = dfh(1:lh,:);
dfm = sum(real(ifft(conj(Fc).*Fha)), 2);dfm = dfm(lm:-1:1,:);

df = [dfh(:); dfm];

% add in huber(D*h_i) \approx ||h_i||_TV term
for i=1:nc
   f = f + lambda_h_TV*huber(D(h(:,i)), huber_d)/huber_d;
   dfh(:,i) = dfh(:,i) + lambda_h_TV*DT(huber_grad(D(h(:,i)), huber_d)/huber_d);
end

end

% "Notes" for computation of Hankel mat-vec
%R = hankel(hc, hr);
%dfh = R*m;
%dfm = R.'*h;

%T = toeplitz(tc, tr); % T == R*fliplr(eye(lm)), P == fliplr(eye(lm))
%dfh = T*m(end:-1:1); % keep track of H*m == (H*P)*P'*m == T*P'*m
%dfm = T.'*h; dfm = dfm(end:-1:1);


% forward and transpose of diff() operator used in 1-D TV norm computation
function [Dx] = D(x)
   Dx = diff(x); % diff along columns
end

function [DTy] = DT(y)
   DTy = zeros(size(y,1)+1, size(y,2));
   DTy(1,:) = -y(1,:);
   DTy(2:end-1,:) = y(1:end-1,:)-y(2:end,:);
   DTy(end,:) = y(end,:);
end

% Huber and friends
% TODO: these aren't vectorized
function [L] = huber(x, d)
   L = sum(0.5*x.^2.*(abs(x) <= d) + d*(abs(x) - 0.5*d).*(abs(x) > d));
end

function [L_grad] = huber_grad(x, d)
   L_grad = x.*(abs(x) <= d) - d*(x < -d) + d*(x > d);
end
