# Adjoint Wavelet
The adjoint (transpose for real wavelets) of wavelet analysis operators appear in some image deblurring problems [1].  As an example, let `R` be a known blurring operator, `W` a wavelet analysis operator, and `b` an image blurred under the action of R.  Since we expect images to have sparse wavlet coefficients, the following problem formulation is reasonable:

```
min_x  ||R*W*x - b||_2^2 + lambda*||x||_1
```

The gradient of the first term is `2*(W^*)*(R^*)*(R*W*x-b)`.  One can use this with (fast) proximal gradient method (e.g. FISTA [1]).  We can compute `R` and `R^*` efficiently in the Fourier domain.  However, wavelet toolboxes do not provide an efficient means to compute `W^*`, as it is not a standard operation.  The standard operations are `W` and `W^+`, the pseudoinverse.

If the wavelet analysis operator is orthogonal, we can reformulate the problem exactly, so we no longer need `W^*` ([2], Proposition 11).  Additionally, we know that `W^*=W^+`, so we could solve the original problem.

For non-orthogonal wavelet operators (e.g. CDF wavelets), the reformulation is not exact.  Our code computes the action of `W^*` in O(N) operations, facilitating the efficient solution of the original problem.

  1. Amir Beck and Marc Teboulle, A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems, SIAM Journal on Imaging Sciences 2 (2009).
  2. Patrick L. Combettes and Jean-Christophe Pesquet, A Douglas–Rachford Splitting Approach to Nonsmooth Convex Variational Signal Recovery, IEEE Journal of Selected Topics in Signal Processing (2007).
