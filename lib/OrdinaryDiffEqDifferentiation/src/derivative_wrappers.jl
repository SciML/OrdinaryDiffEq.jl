const FIRST_AUTODIFF_TGRAD_MESSAGE = """
First call to automatic differentiation for time gradient
failed. This means that the user `f` function is not compatible
with automatic differentiation. Methods to fix this include:

1. Turn off automatic differentiation (e.g. Rosenbrock23() becomes
   Rosenbrock23(autodiff=AutoFiniteDiff())). More details can be found at
   https://docs.sciml.ai/DiffEqDocs/stable/features/performance_overloads/
2. Improving the compatibility of `f` with ForwardDiff.jl automatic
   differentiation (using tools like PreallocationTools.jl). More details
   can be found at https://docs.sciml.ai/DiffEqDocs/stable/basics/faq/#Autodifferentiation-and-Dual-Numbers
3. Defining analytical Jacobians and time gradients. More details can be
   found at https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction

Note 1: this failure occurred inside of the time gradient function. These
time gradients are only required by Rosenbrock methods (`Rosenbrock23`,
`Rodas4`, etc.) and are done by automatic differentiation w.r.t. the
argument `t`. If your function is compatible with automatic differentiation
w.r.t. `u`, i.e. for Jacobian generation, another way to work around this
issue is to switch to a non-Rosenbrock method.

Note 2: turning off automatic differentiation tends to have a very minimal
performance impact (for this use case, because it's forward mode for a
square Jacobian. This is different from optimization gradient scenarios).
However, one should be careful as some methods are more sensitive to
accurate gradients than others. Specifically, Rodas methods like `Rodas4`
and `Rodas5P` require accurate Jacobians in order to have good convergence,
while many other methods like BDF (`QNDF`, `FBDF`), SDIRK (`KenCarp4`),
and Rosenbrock-W (`Rosenbrock23`) do not. Thus if using an algorithm which
is sensitive to autodiff and solving at a low tolerance, please change the
algorithm as well.
"""

struct FirstAutodiffTgradError <: Exception
    e::Any
end

function Base.showerror(io::IO, e::FirstAutodiffTgradError)
    println(io, FIRST_AUTODIFF_TGRAD_MESSAGE)
    Base.showerror(io, e.e)
end

const FIRST_AUTODIFF_JAC_MESSAGE = """
First call to automatic differentiation for the Jacobian
failed. This means that the user `f` function is not compatible
with automatic differentiation. Methods to fix this include:

1. Turn off automatic differentiation (e.g. Rosenbrock23() becomes
   Rosenbrock23(autodiff = AutoFiniteDiff())). More details can befound at
   https://docs.sciml.ai/DiffEqDocs/stable/features/performance_overloads/
2. Improving the compatibility of `f` with ForwardDiff.jl automatic
   differentiation (using tools like PreallocationTools.jl). More details
   can be found at https://docs.sciml.ai/DiffEqDocs/stable/basics/faq/#Autodifferentiation-and-Dual-Numbers
3. Defining analytical Jacobians. More details can be
   found at https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction

Note: turning off automatic differentiation tends to have a very minimal
performance impact (for this use case, because it's forward mode for a
square Jacobian. This is different from optimization gradient scenarios).
However, one should be careful as some methods are more sensitive to
accurate gradients than others. Specifically, Rodas methods like `Rodas4`
and `Rodas5P` require accurate Jacobians in order to have good convergence,
while many other methods like BDF (`QNDF`, `FBDF`), SDIRK (`KenCarp4`),
and Rosenbrock-W (`Rosenbrock23`) do not. Thus if using an algorithm which
is sensitive to autodiff and solving at a low tolerance, please change the
algorithm as well.
"""

struct FirstAutodiffJacError <: Exception
    e::Any
end

function Base.showerror(io::IO, e::FirstAutodiffJacError)
    println(io, FIRST_AUTODIFF_JAC_MESSAGE)
    Base.showerror(io, e.e)
end

function jacobian(f::F, x::AbstractArray{<:Number}, integrator) where {F}
    alg = unwrap_alg(integrator, true)

    # Update stats.nf

    dense = ADTypes.dense_ad(alg_autodiff(alg))

    if dense isa AutoForwardDiff
        sparsity, colorvec = sparsity_colorvec(integrator.f, x)
        maxcolor = maximum(colorvec)
        chunk_size = (get_chunksize(alg) == Val(0) || get_chunksize(alg) == Val(nothing)) ?
            nothing : get_chunksize(alg)
        num_of_chunks = div(
            maxcolor,
            isnothing(chunk_size) ?
                getsize(ForwardDiff.pickchunksize(maxcolor)) : _unwrap_val(chunk_size),
            RoundUp
        )

        integrator.stats.nf += num_of_chunks

    elseif dense isa AutoFiniteDiff
        sparsity, colorvec = sparsity_colorvec(integrator.f, x)
        if dense.fdtype == Val(:forward)
            integrator.stats.nf += maximum(colorvec) + 1
        elseif dense.fdtype == Val(:central)
            integrator.stats.nf += 2 * maximum(colorvec)
        elseif dense.fdtype == Val(:complex)
            integrator.stats.nf += maximum(colorvec)
        end
    else
        integrator.stats.nf += 1
    end

    if dense isa AutoFiniteDiff
        dense = SciMLBase.@set dense.dir = diffdir(integrator)
    end

    # Apply GPU-safe wrapping for AutoForwardDiff when dealing with GPU arrays
    dense = gpu_safe_autodiff(dense, x)

    autodiff_alg = gpu_safe_autodiff(alg_autodiff(alg), x)

    if alg_autodiff(alg) isa AutoSparse
        autodiff_alg = SciMLBase.@set autodiff_alg.dense_ad = dense
    else
        autodiff_alg = dense
    end

    if integrator.iter == 1
        try
            jac = DI.jacobian(f, autodiff_alg, x)
        catch e
            throw(FirstAutodiffJacError(e))
        end
    else
        jac = DI.jacobian(f, autodiff_alg, x)
    end

    return jac
end

# fallback for scalar x, is needed for calc_J to work
function jacobian(f::F, x, integrator) where {F}
    alg = unwrap_alg(integrator, true)

    dense = ADTypes.dense_ad(alg_autodiff(alg))

    if dense isa AutoForwardDiff
        integrator.stats.nf += 1
    elseif dense isa AutoFiniteDiff
        if dense.fdtype == Val(:forward)
            integrator.stats.nf += 2
        elseif dense.fdtype == Val(:central)
            integrator.stats.nf += 2
        elseif dense.fdtype == Val(:complex)
            integrator.stats.nf += 1
        end
    else
        integrator.stats.nf += 1
    end

    if dense isa AutoFiniteDiff
        dense = SciMLBase.@set dense.dir = diffdir(integrator)
    end

    # Apply GPU-safe wrapping for AutoForwardDiff when dealing with GPU arrays
    dense = gpu_safe_autodiff(dense, x)

    autodiff_alg = gpu_safe_autodiff(alg_autodiff(alg), x)

    if autodiff_alg isa AutoSparse
        autodiff_alg = SciMLBase.@set autodiff_alg.dense_ad = dense
    else
        autodiff_alg = dense
    end

    if integrator.iter == 1
        try
            jac = DI.derivative(f, autodiff_alg, x)
        catch e
            throw(FirstAutodiffJacError(e))
        end
    else
        jac = DI.derivative(f, autodiff_alg, x)
    end

    return jac
end

function jacobian!(
        J::AbstractMatrix{<:Number}, f::F, x::AbstractArray{<:Number},
        fx::AbstractArray{<:Number}, integrator::SciMLBase.DEIntegrator,
        jac_config
    ) where {F}
    # Handle empty state vector - nothing to compute
    if isempty(x)
        return nothing
    end

    alg = unwrap_alg(integrator, true)

    dense = ADTypes.dense_ad(alg_autodiff(alg))

    if dense isa AutoForwardDiff
        if alg_autodiff(alg) isa AutoSparse
            integrator.stats.nf += maximum(ncolors(jac_config[1]))
        else
            sparsity, colorvec = sparsity_colorvec(integrator.f, x)
            maxcolor = maximum(colorvec)
            chunk_size = (
                    get_chunksize(alg) == Val(0) ||
                    get_chunksize(alg) == Val(nothing)
                ) ? nothing : get_chunksize(alg)
            num_of_chunks = chunk_size === nothing ?
                Int(
                    ceil(
                        maxcolor /
                        getsize(ForwardDiff.pickchunksize(maxcolor))
                    )
                ) :
                Int(ceil(maxcolor / _unwrap_val(chunk_size)))

            integrator.stats.nf += num_of_chunks
        end

    elseif dense isa AutoFiniteDiff
        sparsity, colorvec = sparsity_colorvec(integrator.f, x)
        if dense.fdtype == Val(:forward)
            integrator.stats.nf += maximum(colorvec) + 1
        elseif dense.fdtype == Val(:central)
            integrator.stats.nf += 2 * maximum(colorvec)
        elseif dense.fdtype == Val(:complex)
            integrator.stats.nf += maximum(colorvec)
        end
    else
        integrator.stats.nf += 1
    end

    if dense isa AutoFiniteDiff
        config = diffdir(integrator) > 0 ? jac_config[1] : jac_config[2]
    else
        config = jac_config[1]
    end

    if integrator.iter == 1
        try
            DI.jacobian!(f, fx, J, config, gpu_safe_autodiff(alg_autodiff(alg), x), x)
        catch e
            throw(FirstAutodiffJacError(e))
        end
    else
        DI.jacobian!(f, fx, J, config, gpu_safe_autodiff(alg_autodiff(alg), x), x)
    end

    return nothing
end

function build_jac_config(
        alg, f::F1, uf::F2, du1, uprev,
        u, tmp, du2
    ) where {F1, F2}
    haslinsolve = hasfield(typeof(alg), :linsolve)

    if !SciMLBase.has_jac(f) &&
            (!SciMLBase.has_Wfact_t(f)) &&
            (
            (
                concrete_jac(alg) === nothing && (
                    !haslinsolve || (
                        haslinsolve &&
                            (alg.linsolve === nothing || LinearSolve.needs_concrete_A(alg.linsolve))
                    )
                )
            ) ||
                (concrete_jac(alg) !== nothing && concrete_jac(alg))
        )
        jac_prototype = f.jac_prototype

        if is_sparse_csc(jac_prototype)
            if f.mass_matrix isa UniformScaling
                idxs = diagind(jac_prototype)
                @. @view(jac_prototype[idxs]) = 1
            else
                idxs = findall(!iszero, f.mass_matrix)
                @. @view(jac_prototype[idxs]) = @view(f.mass_matrix[idxs])
            end
        end

        autodiff_alg = gpu_safe_autodiff(alg_autodiff(alg), u)
        dense = autodiff_alg isa AutoSparse ? ADTypes.dense_ad(autodiff_alg) : autodiff_alg

        if dense isa AutoFiniteDiff
            dir_forward = @set dense.dir = 1
            dir_reverse = @set dense.dir = -1

            if autodiff_alg isa AutoSparse
                autodiff_alg_forward = @set autodiff_alg.dense_ad = dir_forward
                autodiff_alg_reverse = @set autodiff_alg.dense_ad = dir_reverse
            else
                autodiff_alg_forward = dir_forward
                autodiff_alg_reverse = dir_reverse
            end

            jac_config_forward = DI.prepare_jacobian(
                uf, du1, autodiff_alg_forward, u, strict = Val(false)
            )
            jac_config_reverse = DI.prepare_jacobian(
                uf, du1, autodiff_alg_reverse, u, strict = Val(false)
            )

            jac_config = (jac_config_forward, jac_config_reverse)
        else
            jac_config1 = DI.prepare_jacobian(uf, du1, autodiff_alg, u, strict = Val(false))
            jac_config = (jac_config1, jac_config1)
        end

    else
        jac_config = (nothing, nothing)
    end

    return jac_config
end

function get_chunksize(
        jac_config::ForwardDiff.JacobianConfig{
            T,
            V,
            N,
            D,
        }
    ) where {
        T, V, N, D,
    }
    return Val(N)
end # don't degrade compile time information to runtime information

function resize_jac_config!(cache, integrator)
    if !isnothing(cache.jac_config) && !isnothing(cache.jac_config[1])
        uf = cache.uf
        uf = SciMLBase.@set uf.f = SciMLBase.unwrapped_f(uf.f)

        # for correct FiniteDiff dirs
        autodiff_alg = alg_autodiff(integrator.alg)
        if autodiff_alg isa AutoFiniteDiff
            ad_right = SciMLBase.@set autodiff_alg.dir = 1
            ad_left = SciMLBase.@set autodiff_alg.dir = -1
        else
            ad_right = autodiff_alg
            ad_left = autodiff_alg
        end

        cache.jac_config = (
            [
                DI.prepare!_jacobian(
                        uf, cache.du1, config, ad, integrator.u
                    )
                    for (ad, config) in zip(
                        (ad_right, ad_left), cache.jac_config
                    )
            ]...,
        )
    end
    return cache.jac_config
end

function resize_grad_config!(cache, integrator)
    if !isnothing(cache.grad_config) && !isnothing(cache.grad_config[1])

        # for correct FiniteDiff dirs
        autodiff_alg = ADTypes.dense_ad(alg_autodiff(integrator.alg))
        if autodiff_alg isa AutoFiniteDiff
            ad_right = SciMLBase.@set autodiff_alg.dir = 1
            ad_left = SciMLBase.@set autodiff_alg.dir = -1
        else
            ad_right = autodiff_alg
            ad_left = autodiff_alg
        end

        cache.grad_config = (
            [
                DI.prepare!_derivative(
                        cache.tf, cache.du1, config, ad, integrator.t
                    )
                    for (ad, config) in zip(
                        (ad_right, ad_left), cache.grad_config
                    )
            ]...,
        )
    end
    return cache.grad_config
end

"""
    gpu_safe_autodiff(backend, u)

Automatically wrap AutoForwardDiff with AutoForwardFromPrimitive for GPU arrays
that don't support fast scalar indexing (e.g., GPU arrays).
"""
function gpu_safe_autodiff(backend::AutoForwardDiff, u)
    if ArrayInterface.fast_scalar_indexing(u)
        # CPU arrays with fast scalar indexing - use original backend
        return backend
    else
        # GPU arrays or arrays without fast scalar indexing - use primitive wrapper
        return DI.AutoForwardFromPrimitive(backend)
    end
end

# Fallback for other AD backends
gpu_safe_autodiff(backend, u) = backend

function build_grad_config(alg, f::F1, tf::F2, du1, t) where {F1, F2}
    if !SciMLBase.has_tgrad(f)
        ad = ADTypes.dense_ad(alg_autodiff(alg))

        # Apply GPU-safe wrapping for AutoForwardDiff when dealing with GPU arrays
        ad = gpu_safe_autodiff(ad, du1)

        if ad isa AutoFiniteDiff
            dir_true = @set ad.dir = 1
            dir_false = @set ad.dir = -1

            grad_config_true = DI.prepare_derivative(tf, du1, dir_true, t)
            grad_config_false = DI.prepare_derivative(tf, du1, dir_false, t)

            grad_config = (grad_config_true, grad_config_false)
        elseif ad isa AutoForwardDiff
            grad_config1 = DI.prepare_derivative(tf, du1, ad, convert(eltype(du1), t))
            grad_config = (grad_config1, grad_config1)
        else
            grad_config1 = DI.prepare_derivative(tf, du1, ad, t)
            grad_config = (grad_config1, grad_config1)
        end
        return grad_config
    else
        return (nothing, nothing)
    end
end

function sparsity_colorvec(f::F, x) where {F}
    sparsity = f.sparsity

    if is_sparse_csc(sparsity)
        if f.mass_matrix isa UniformScaling
            idxs = diagind(sparsity)
            @. @view(sparsity[idxs]) = 1
        else
            idxs = findall(!iszero, f.mass_matrix)
            @. @view(sparsity[idxs]) = @view(f.mass_matrix[idxs])
        end
    end

    col_alg = GreedyColoringAlgorithm()
    col_prob = ColoringProblem()
    colorvec = SciMLBase.has_colorvec(f) ? f.colorvec :
        (
            isnothing(sparsity) ? (1:length(x)) :
            column_colors(coloring(sparsity, col_prob, col_alg))
        )
    return sparsity, colorvec
end
