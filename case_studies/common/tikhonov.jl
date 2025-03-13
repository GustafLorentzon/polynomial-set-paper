# Create an approximate Tikhonov solution to
# the linear system Jx=f with regularization parameter lambda

function tikhonov_step(J, f, lambda)
    # Compute the SVD of J
    U, S, V = svd(J)

    # Compute the Tikhonov filter factors
    S_reg = Diagonal(S ./ (S.^2 .+ lambda^2))

    # Compute the regularized solution
    x_lambda = V * (S_reg * (U' * f))

    return x_lambda
end
