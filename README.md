# Adjoint Wavelet
The adjoint (transpose for real wavelets) of wavelet analysis operators appear in some image deblurring problems [1].  As an example, let `R` be a known blurring operator, `W` a wavelet analysis operator, and `b` an image blurred under the action of R.  Since we expect images to have sparse wavlet coefficients, the following problem formulation is reasonable:

```
min_x  ||R*W*x - b||_2^2 + lambda*||x||_1
```

The gradient of the first term is `2*(W^*)*(R^*)*(R*W*x-b)`.  One can use this with (fast) proximal gradient method (e.g. FISTA [1]).  We can compute `R` and `R^*` efficiently in the Fourier domain.  However, wavelet toolboxes do not provide an efficient means to compute `W^*`, as it is not a standard operation.  The standard operations are `W` and `W^+`, the pseudoinverse.

If the wavelet analysis operator is orthogonal, we can reformulate the problem exactly, so we no longer need `W^*` ([2], Proposition 11).  Additionally, we know that `W^*=W^+`, so we could solve the original problem.

For non-orthogonal wavelet operators (e.g. CDF wavelets), the reformulation is not exact.  Our code computes the action of `W^*` in O(N) operations, facilitating the efficient solution of the original problem.

## MATLAB
Most of the code is in MATLAB.  I've also put [Beck's FISTA code](https://web.iem.technion.ac.il/images/user-files/becka/papers/wavelet_FISTA.zip) and the [HNO image deblurring code](http://www.imm.dtu.dk/~pcha/HNO/) in the repo for some of our experiments.

## Julia
There is a bit of experimental code in the julia directory.  We tried to see how fast automatic differentiation could compute the gradient.


# Blind Channel Estimation - Unconstrained Formulation
Suppose a single unknown source signal is sent over a few channels with unknown impulse responses.  We observe and record the output of each channel and wish to recover the source signal and channel impulse responses.  We pose an unconstrained non-convex problem and solve it with Mark Schmidt's L1General.  To use the routines in L1General, we must provide it with a gradient subroutine.  The gradients of the differentiable terms in our formulation involve an interesting adjoint, which we show how to compute efficiently.


# IEEE SPM
We're working on a ``lecture note'' for IEEE's Signal Processing Magazine that describes a framework for computing the adjoint of fast discrete wavelet transform and also presents the BCE problem and the interesting adjoint.


# SIAM Front Range Applied Math Student Conference
I gave a quick talk at the 2016 SIAM FRAMSC in Denver, CO.  There are a bunch of slides, many of which I didn't present!


  1. Amir Beck and Marc Teboulle, *A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems*, SIAM Journal on Imaging Sciences 2 (2009).
  2. Patrick L. Combettes and Jean-Christophe Pesquet, *A Douglas–Rachford Splitting Approach to Nonsmooth Convex Variational Signal Recovery*, IEEE Journal of Selected Topics in Signal Processing (2007).
