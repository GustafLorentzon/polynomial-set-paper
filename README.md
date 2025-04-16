# polynomial-set-paper
Code associated with the paper “The polynomial set associated with a fixed number of matrix-matrix multiplications” (Elias Jarlebring and Gustaf Lorentzon, 2025):
[https://arxiv.org/abs/2504.01500](https://arxiv.org/abs/2504.01500)

./common: Contains shared files used when generating evaluation scheme coefficients.

./data: Stores all precomputed evaluation scheme coefficients.

./five_mult Contains files for generating evaluation schemes that compute degree-20 matrix polynomials using only 5 matrix-matrix multiplications.

./six_mult Contains files for generating schemes that compute truncated Taylor polynomials of the exponential function — achieving degree-30 (real coefficients) and degree-32 (complex coefficients) — using only 6 matrix-matrix multiplications.

./seven_mult Contains files for generating schemes that compute degree-42 polynomials (complex coefficients) using only 7 matrix-matrix multiplications.

dimension_verificiation.jl Polts the singular values of the Jacobian in a random point in two different precisions, allowing us to count the number of non-zero singular values, even when the jacobian is very ill-conditioned.

reduction_figure.jl Plots the degrees of polynomials corresponding to different reduction structures in the evaluation scheme, see Figure 4.1 in [https://arxiv.org/abs/2504.01500](https://arxiv.org/abs/2504.01500).
