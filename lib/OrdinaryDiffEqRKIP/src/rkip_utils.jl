# cannot use @bb on axpy due to https://github.com/SciML/MaybeInplace.jl/issues/14
@inline axpy_mip(α, x::uType, y::uType, ::Val{true}) where {uType} = axpy!(α, x, y) #BLAS call
@inline axpy_mip(α, x::uType, y::uType, ::Val{false}) where {uType} = y .+ α .* x

"""
Maybe-In Place "mip"
Safe function for computing the matrix exponential vector product for both mutable and immutable type.
Same principle of operation as `_safe_matvec_prod` for mutability/in-place handling.
"""
@inline expmv_rkip_mip(cache::RKIPCache{expOpType, cacheType, tType, opType, uType, iip}, v::uType, h::tType, p, t) where {expOpType, cacheType, tType, opType, uType, iip} = expmv_rkip_mip(
    cache, v, h, lastindex(cache.c_mapping), p, t) # If i is not precised, this is the final step in the RKIP -> h = dt

@inline function expmv_rkip_mip(
        cache::RKIPCache{expOpType, cacheType, tType, opType, uType, iip}, v::uType,
        h::tType, stage_index::Int, p, t) where {
        expOpType, cacheType, tType, opType, uType, iip}
    if !(h ≈ tType(0.0))
        c = cache.c_mapping[stage_index]
        exp_cache = cache.exp_cache
        tmp = cache.tmp
        return _expmv_rkip_mip(exp_cache, tmp, v, c, h > tType(0.0), iip, p, t)
    end
    return v
end

@inline function _expmv_rkip_mip(
        cache::ExpCacheNoLdiv{expOpType}, tmp::vType, v::vType, stage_index::Integer,
        positive, iip, p, t) where {expOpType <: AbstractSciMLOperator, vType}
    op = get_op_for_this_step(cache, positive, stage_index)
    v = matvec_prod_mip(tmp, op, v, iip, p, t)
    return v
end

@inline function _expmv_rkip_mip(
        cache::ExpCache{expOpType}, tmp::vType, v::vType, stage_index::Integer,
        positive, ::Val{true}, p, t) where {expOpType <: AbstractSciMLOperator, vType}
    op = get_op_for_this_step(cache, stage_index)
    if positive
        matvec_prod_mip(tmp, op, v, Val(true), p, t)
    else
        ldiv!(tmp, op, v) #ldiv can only be used for in place call
        copyto!(v, tmp)
    end
    return v
end

"""
Maybe-In Place "mip"
Safe function for matrix vector product for both mutable and immutable type.
For mutable type, overwrite `v` with `mat*v` and return `v`. Otherwise return  `v`
Use the dispatch between ::Val{true} (in place) and ::Val{false} immutable to decide
"""
@inline function matvec_prod_mip(
        tmp::V, mat::M, v::V, ::Val{true}, p, t) where {V, M <: AbstractSciMLOperator}
    mat(tmp, v, v, p, t)
    copyto!(v, tmp)
    return v
end

@inline function matvec_prod_mip(
        _::V, mat::M, v::V, ::Val{false}, p, t) where {V, M <: AbstractSciMLOperator}
    v = mat(v, v, p, t)
    return v
end
