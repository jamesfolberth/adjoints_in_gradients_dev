
using PyPlot
using Wavelets
using ForwardDiff

srand(0)

#wt_g = wavelet(WT.haar, WT.Filter)
wt_g = wavelet(WT.cdf97, WT.Lifting)

function check_adjoint()

   N = 8
   L = 1
   I = eye(N)
   W = zeros(N,N)
   Wt = zeros(N,N)

   for i in 1:N
      W[:,i] = dwt(I[:,i], wt_g)
      Wt[:,i] = idwt(I[:,i], wt_g)
   end
   
   dump(W)
   dump(Wt)

   println(norm(W.'-Wt,Inf))

end

check_adjoint()

