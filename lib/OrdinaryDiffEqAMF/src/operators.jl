"""
    AMFOperator(jac)
    AMFOperator(split_ops)
    AMFOperator(; factors)

Constructs the operator used as the `W_prototype` for Rosenbrock-W methods.

- `AMFOperator(jac)` builds the exact Rosenbrock operator `-(I - gamma * J) / gamma`.
- `AMFOperator(split_ops)` builds the AMF approximation from an ordered collection of split
  Jacobian operators, `-prod(I - gamma * J_i) / gamma`.
- `AMFOperator(; factors)` assembles `-prod(factors) / gamma` from pre-built AMF factors.
"""
function AMFOperator(jac::SciMLOperators.AbstractSciMLOperator)
    N = size(jac, 1)
    size(jac, 1) == size(jac, 2) || throw(ArgumentError("Jacobian operators must be square"))

    gamma_op = ScalarOperator(
        1.0;
        update_func = (old_op, u, p, t; gamma = 1.0) -> gamma,
        accepted_kwargs = Val((:gamma,)),
    )
    transform_op = ScalarOperator(
        0.0;
        update_func = (old_op, u, p, t; gamma = 1.0) -> inv(gamma),
        accepted_kwargs = Val((:gamma,)),
    )

    return -(IdentityOperator(N) - gamma_op * jac) * transform_op
end

function AMFOperator(split_ops::Union{Tuple, AbstractVector})
    split_tuple = _normalize_split(split_ops)
    _validate_split_sizes(split_tuple)
    factors = map(_amf_factor_operator, split_tuple)
    return _assemble_amf_operator(factors)
end

function AMFOperator(; factors)
    factor_tuple = _normalize_split(factors)
    _validate_split_sizes(factor_tuple)
    return _assemble_amf_operator(factor_tuple)
end

function _normalize_split(split_ops::Union{Tuple, AbstractVector})
    split_tuple = Tuple(split_ops)
    isempty(split_tuple) && throw(ArgumentError("split must contain at least one operator"))
    return split_tuple
end

function _validate_split_sizes(split_ops::Tuple)
    first_op = first(split_ops)
    N = size(first_op, 1)
    size(first_op, 1) == size(first_op, 2) ||
        throw(ArgumentError("split operators must be square"))

    for op in Base.tail(split_ops)
        size(op, 1) == N && size(op, 2) == N ||
            throw(ArgumentError("all split operators must have the same square size"))
    end

    return N
end

function _amf_factor_operator(J_op::SciMLOperators.AbstractSciMLOperator)
    N = size(J_op, 1)
    size(J_op, 1) == size(J_op, 2) || throw(ArgumentError("Jacobian operators must be square"))

    gamma_op = ScalarOperator(
        1.0;
        update_func = (old_op, u, p, t; gamma = 1.0) -> gamma,
        accepted_kwargs = Val((:gamma,)),
    )

    return IdentityOperator(N) - gamma_op * J_op
end

function _assemble_amf_operator(factors)
    transform_op = ScalarOperator(
        0.0;
        update_func = (old_op, u, p, t; gamma = 1.0) -> inv(gamma),
        accepted_kwargs = Val((:gamma,)),
    )

    factor_product = reduce(*, factors)
    return -(factor_product) * transform_op
end

function build_amf_function(
        f!;
        jac,
        split = nothing,
        amf_factors = nothing,
        jac_cache = nothing,
        w_cache = nothing,
        sparsity = nothing,
    )

    jac isa SciMLOperators.AbstractSciMLOperator ||
        throw(ArgumentError("jac must be an AbstractSciMLOperator"))

    N = size(jac, 1)
    size(jac, 1) == size(jac, 2) || throw(ArgumentError("jac must be square"))

    isnothing(jac_cache) && (jac_cache = zeros(N))
    isnothing(w_cache) && (w_cache = zeros(N))

    J_op = cache_operator(jac, jac_cache)

    if !isnothing(split) && !isnothing(amf_factors)
        split_tuple = _normalize_split(split)
        factor_tuple = _normalize_split(amf_factors)
        length(split_tuple) == length(factor_tuple) ||
            throw(ArgumentError("split and amf_factors must contain the same number of operators"))
    end

    W_base = if !isnothing(amf_factors)
        factor_tuple = _normalize_split(amf_factors)
        _validate_split_sizes(factor_tuple) == N ||
            throw(ArgumentError("amf_factors must have the same size as jac"))
        AMFOperator(; factors = factor_tuple)
    elseif isnothing(split)
        AMFOperator(J_op)
    else
        split_tuple = _normalize_split(split)
        _validate_split_sizes(split_tuple) == N ||
            throw(ArgumentError("split operators must have the same size as jac"))
        AMFOperator(split_tuple)
    end
    W_op = cache_operator(W_base, w_cache)

    actual_sparsity = isnothing(sparsity) ? convert(AbstractMatrix, J_op) : sparsity
    return SciMLBase.ODEFunction(f!; jac_prototype = J_op, W_prototype = W_op, sparsity = actual_sparsity)
end
