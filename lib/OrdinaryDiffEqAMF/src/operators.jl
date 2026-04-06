"""
    AMFOperator(J1_op, J2_op)

Constructs an Approximate Matrix Factorization (AMF) operator for implicit time integration.

# Theory

In implicit methods, we need to solve systems of the form:
```
(I - γJ) * Δu = rhs
```

where `J` is the Jacobian. The AMF approximation splits the Jacobian as:
```
J = J1 + J2
```

and approximates:
```
(I - γJ)^(-1) ≈ (I - γJ2)^(-1) * (I - γJ1)^(-1)
```

This is particularly efficient when J1 and J2 are easily invertible (e.g., upper/lower triangular,
diagonal, etc.).

# Arguments
- `J1_op::AbstractSciMLOperator`: First Jacobian split operator (represents J1 only)
- `J2_op::AbstractSciMLOperator`: Second Jacobian split operator (represents J2 only)

# Returns
An operator W representing the AMF approximation: `W = -(I - γJ1)(I - γJ2) / γ`
"""
function AMFOperator(
        J1_op::SciMLOperators.AbstractSciMLOperator,
        J2_op::SciMLOperators.AbstractSciMLOperator,
    )
    N = size(J1_op, 1)
    @assert size(J1_op) == size(J2_op) "J1 and J2 must have the same size"

    I_N = Diagonal(ones(N))

    W1_prototype = convert(AbstractMatrix, I_N - J1_op)

    W1_op = MatrixOperator(
        W1_prototype;
        update_func! = (M, u, p, t; gamma = 1.0) -> begin
            update_coefficients!(J1_op, u, p, t)
            J1_mat = J1_op.A
            copyto!(M, J1_mat)
            lmul!(gamma, M)
            @. M = I_N - M
            nothing
        end,
        accepted_kwargs = Val((:gamma,)),
    )

    W2_prototype = convert(AbstractMatrix, I_N - J2_op)

    W2_op = MatrixOperator(
        W2_prototype;
        update_func! = (M, u, p, t; gamma = 1.0) -> begin
            update_coefficients!(J2_op, u, p, t)
            J2_mat = J2_op.A
            copyto!(M, J2_mat)
            lmul!(gamma, M)
            @. M = I_N - M
            nothing
        end,
        accepted_kwargs = Val((:gamma,)),
    )

    transform_op = ScalarOperator(
        0.0;
        update_func = (old_op, u, p, t; gamma = 1.0) -> inv(gamma),
        accepted_kwargs = Val((:gamma,)),
    )

    W = -(W1_op * W2_op) * transform_op

    return W
end

"""
    AMFOperator(J1_prototype, J2_prototype, fjac_J1!, fjac_J2!)

Constructs an AMF operator using direct Jacobian update functions.

This version calls the Jacobian update functions directly on the W1 and W2 operator matrices,
avoiding indirection through operator composition.

# Arguments
- `J1_prototype::AbstractMatrix`: Prototype matrix for J1 (e.g., `UpperTriangular(zeros(N,N))`)
- `J2_prototype::AbstractMatrix`: Prototype matrix for J2 (e.g., `LowerTriangular(zeros(N,N))`)
- `fjac_J1!::Function`: Function `(M, u, p, t) -> nothing` that fills M with J1 values
- `fjac_J2!::Function`: Function `(M, u, p, t) -> nothing` that fills M with J2 values

# Returns
An operator W representing: `W = -(I - γJ1)(I - γJ2) / γ`
"""
function AMFOperator(
        J1_prototype::AbstractMatrix,
        J2_prototype::AbstractMatrix,
        fjac_J1!::Function,
        fjac_J2!::Function,
    )
    N = size(J1_prototype, 1)
    @assert size(J1_prototype) == size(J2_prototype) "J1 and J2 must have the same size"
    @assert N == size(J1_prototype, 2) "J1 and J2 must be square"

    I_N = Diagonal(ones(N))

    W1_prototype = convert(AbstractMatrix, I_N - J1_prototype)

    W1_op = MatrixOperator(
        W1_prototype;
        update_func! = (M, u, p, t; gamma = 1.0) -> begin
            fjac_J1!(M, u, p, t)
            lmul!(gamma, M)
            @. M = I_N - M
            nothing
        end,
        accepted_kwargs = Val((:gamma,)),
    )

    W2_prototype = convert(AbstractMatrix, I_N - J2_prototype)

    W2_op = MatrixOperator(
        W2_prototype;
        update_func! = (M, u, p, t; gamma = 1.0) -> begin
            fjac_J2!(M, u, p, t)
            lmul!(gamma, M)
            @. M = I_N - M
            nothing
        end,
        accepted_kwargs = Val((:gamma,)),
    )

    transform_op = ScalarOperator(
        0.0;
        update_func = (old_op, u, p, t; gamma = 1.0) -> inv(gamma),
        accepted_kwargs = Val((:gamma,)),
    )

    W = -(W1_op * W2_op) * transform_op

    return W
end

function split_jacobian_operator(
        n::Integer,
        jac_upper,
        jac_lower;
        u_cache = zeros(n^2),
    )
    J1_op = MatrixOperator(UpperTriangular(zeros(n, n)); update_func! = jac_upper)
    J2_op = MatrixOperator(LowerTriangular(zeros(n, n)); update_func! = jac_lower)
    return cache_operator(J1_op + J2_op, u_cache)
end

function build_amf_function(
        f!;
        n::Integer,
        jac_upper,
        jac_lower,
        jac_cache = zeros(n^2),
        w_cache = zeros(n),
    )
    J_op = split_jacobian_operator(n, jac_upper, jac_lower; u_cache = jac_cache)
    W_op = AMFOperator(
        UpperTriangular(zeros(n, n)),
        LowerTriangular(zeros(n, n)),
        jac_upper,
        jac_lower,
    )
    W_op = cache_operator(W_op, w_cache)
    return SciMLBase.ODEFunction(f!; jac_prototype = J_op, W_prototype = W_op, sparsity = convert(AbstractMatrix, J_op))
end
