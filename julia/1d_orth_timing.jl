
using PyPlot
using Wavelets
using ForwardDiff

srand(0)

wt_g = wavelet(WT.haar, WT.Filter)
#wt_g = wavelet(WT.haar, WT.Lifting)

function test_me{T}(x::Array{T,1})
   return sum(sin, x) + prod(tan, x) * sum(sqrt, x)
   #return sum(sin(x)) + prod(tan(x))*sum(sqrt(x))
end


"""
0.5*\|W*x\|_2^2

See Subsection 1.6 of Beck's book chapter on gradient-based algorithms for signal recovery.
The problem posed there treats W as an inverse wavelet transform, so x contains the wavelet coefficients.
This shouldn't be important for orthogonal wavelets.
"""
function f{T}(x::Array{T,1})
   xt = idwt(x, wt_g, 3)
   return 0.5*sumabs2(xt) 
end


"""
W^T*W*x
"""
function Df{T}(x::Array{T,1})
   Wx = idwt(x, wt_g, 3)
   WtWx = dwt(Wx, wt_g, 3)
   return WtWx
end


"""
A^T*A*x

This is just used to time dense mat-vec
"""
function mat_vec{T}(x::Array{T,1}, A::Array{T,2})
   Ax = A*x
   AtAx = A.'*Ax
   return AtAx
end

function example_1()
   
   x = randn(2^10);
   
   Df_AD = gradient(f)
  
   println("Df:")
   @time Df(x)
   println("Df_AD:")
   @time Df_AD(x)
   println("norm(Df(x) - Df_AD(x)) = $(norm(Df(x) - Df_AD(x)))")

end


function runtime_scaling(;dyadic=false)
   
   n_samples = 3
   
   if dyadic
      Nv = 2.^(3:1:14)
   else
      Nv = Array{Int64}(2^3*round(collect(logspace(0,3,10))))
   end
   
   Df_AD = gradient(f)
   
   Df_times = zeros(size(Nv))
   Df_AD_times = zeros(size(Nv))
   mat_vec_times = zeros(size(Nv))

   for (i,N) in enumerate(Nv)
      println(N)
      x = randn(N)
      
      # force compilation
      Df1 = Df(x)
      Df2 = Df_AD(x)
      #A = randn(length(x), length(x))
      A = idwt(eye(length(x)), wt_g, 3)
      Df3 = mat_vec(x,A)
      println("norm(Df(x) - Df_AD(x)) = $(norm(Df1 - Df2))")
      println("norm(Df(x) - mat_vec(x)) = $(norm(Df1 - Df3))")
   
      for s in 1:n_samples
         Df_timed = @timed Df(x)
         Df_AD_timed = @timed Df_AD(x)
         mat_vec_timed = @timed mat_vec(x,A)
      
         Df_times[i] += Df_timed[2]/n_samples
         Df_AD_times[i] += Df_AD_timed[2]/n_samples
         mat_vec_times[i] += mat_vec_timed[2]/n_samples
      end
   end
   
   figure("1")
   clf()
   hold(true)
   loglog(Nv, Df_times, "b-", Nv, Df_AD_times, "g-", Nv, mat_vec_times, "r-", Nv, Df_times, "bo", Nv, Df_AD_times, "go", Nv, mat_vec_times, "ro")
   hold(false)
   xlabel("N")
   ylabel("Time (seconds)")
   legend(["Df", "Df_AD", "mat_vec"], loc="upper left")
   
   if dyadic
      title("Wavelet Transform: Dyadic N")
   else
      title("Wavelet Transform: Non-dyadic N")
   end
end

#example_1()
runtime_scaling(dyadic=false)

