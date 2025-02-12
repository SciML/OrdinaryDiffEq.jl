const FIRST_AUTODIFF_TGRAD_MESSAGE = """
                               First call to automatic differentiation for time gradient
                               failed. This means that the user `f` function is not compatible
                               with automatic differentiation. Methods to fix this include:

                               1. Turn off automatic differentiation (e.g. Rosenbrock23() becomes
                                  Rosenbrock23(autodiff=false)). More details can be found at
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
                                  Rosenbrock23(autodiff=false)). More details can befound at
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

function derivative!(df::AbstractArray{<:Number}, f,
        x::Union{Number, AbstractArray{<:Number}}, fx::AbstractArray{<:Number},
        integrator, grad_config)
    alg = unwrap_alg(integrator, true)
    tmp = length(x) # We calculate derivative for all elements in gradient
    autodiff_alg = alg_autodiff(alg)
    if autodiff_alg isa AutoForwardDiff
        T = if standardtag(alg)
            typeof(ForwardDiff.Tag(OrdinaryDiffEqTag(), eltype(df)))
        else
            typeof(ForwardDiff.Tag(f, eltype(df)))
        end

        xdual = Dual{T, eltype(df), 1}(convert(eltype(df), x),
            ForwardDiff.Partials((one(eltype(df)),)))

        if integrator.iter == 1
            try
                f(grad_config, xdual)
            catch e
                throw(FirstAutodiffTgradError(e))
            end
        else
            f(grad_config, xdual)
        end

        df .= first.(ForwardDiff.partials.(grad_config))
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    elseif autodiff_alg isa AutoFiniteDiff
        FiniteDiff.finite_difference_gradient!(df, f, x, grad_config,
            dir = diffdir(integrator))
        fdtype = alg_difftype(alg)
        if fdtype == Val{:forward} || fdtype == Val{:central}
            tmp *= 2
            if eltype(df) <: Complex
                tmp *= 2
            end
        end
        integrator.stats.nf += tmp
    else
        error("$alg_autodiff not yet supported in derivative! function")
    end
    nothing
end

function derivative(f, x::Union{Number, AbstractArray{<:Number}},
        integrator)
    local d
    tmp = length(x) # We calculate derivative for all elements in gradient
    alg = unwrap_alg(integrator, true)
    if alg_autodiff(alg) isa AutoForwardDiff
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        if integrator.iter == 1
            try
                d = ForwardDiff.derivative(f, x)
            catch e
                throw(FirstAutodiffTgradError(e))
            end
        else
            d = ForwardDiff.derivative(f, x)
        end
    elseif alg_autodiff(alg) isa AutoFiniteDiff
        d = FiniteDiff.finite_difference_derivative(f, x, alg_difftype(alg),
            dir = diffdir(integrator))
        if alg_difftype(alg) === Val{:central} || alg_difftype(alg) === Val{:forward}
            tmp *= 2
        end
        integrator.stats.nf += tmp
        d
    else
        error("$alg_autodiff not yet supported in derivative function")
    end
end

jacobian_autodiff(f, x, odefun, alg) = (ForwardDiff.derivative(f, x), 1, alg)
function jacobian_autodiff(f, x::AbstractArray, odefun, alg)
    jac_prototype = odefun.jac_prototype
    sparsity, colorvec = sparsity_colorvec(odefun, x)
    maxcolor = maximum(colorvec)
    chunk_size = get_chunksize(alg) === Val(0) ? nothing : get_chunksize(alg)
    num_of_chunks = chunk_size === nothing ?
                    Int(ceil(maxcolor / getsize(ForwardDiff.pickchunksize(maxcolor)))) :
                    Int(ceil(maxcolor / _unwrap_val(chunk_size)))
    (
        forwarddiff_color_jacobian(f, x, colorvec = colorvec, sparsity = sparsity,
            jac_prototype = jac_prototype, chunksize = chunk_size),
        num_of_chunks)
end

function _nfcount(N, ::Type{diff_type}) where {diff_type}
    if diff_type === Val{:complex}
        tmp = N
    elseif diff_type === Val{:forward}
        tmp = N + 1
    else
        tmp = 2N
    end
    tmp
end

function jacobian_finitediff(f, x, ::Type{diff_type}, dir, colorvec, sparsity,
        jac_prototype) where {diff_type}
    (FiniteDiff.finite_difference_derivative(f, x, diff_type, eltype(x), dir = dir), 2)
end
function jacobian_finitediff(f, x::AbstractArray, ::Type{diff_type}, dir, colorvec,
        sparsity, jac_prototype) where {diff_type}
    f_in = diff_type === Val{:forward} ? f(x) : similar(x)
    ret_eltype = eltype(f_in)
    J = FiniteDiff.finite_difference_jacobian(f, x, diff_type, ret_eltype, f_in,
        dir = dir, colorvec = colorvec,
        sparsity = sparsity,
        jac_prototype = jac_prototype)
    return J, _nfcount(maximum(colorvec), diff_type)
end
function jacobian(f, x, integrator)
    alg = unwrap_alg(integrator, true)
    local tmp
    if alg_autodiff(alg) isa AutoForwardDiff
        if integrator.iter == 1
            try
                J, tmp = jacobian_autodiff(f, x, integrator.f, alg)
            catch e
                throw(FirstAutodiffJacError(e))
            end
        else
            J, tmp = jacobian_autodiff(f, x, integrator.f, alg)
        end
    elseif alg_autodiff(alg) isa AutoFiniteDiff
        jac_prototype = integrator.f.jac_prototype
        sparsity, colorvec = sparsity_colorvec(integrator.f, x)
        dir = diffdir(integrator)
        J, tmp = jacobian_finitediff(f, x, alg_difftype(alg), dir, colorvec, sparsity,
            jac_prototype)
    else
        bleh
    end
    integrator.stats.nf += tmp
    J
end

function jacobian_finitediff_forward!(J, f, x, jac_config, forwardcache, integrator)
    (FiniteDiff.finite_difference_jacobian!(J, f, x, jac_config, forwardcache,
        dir = diffdir(integrator));
    maximum(jac_config.colorvec))
end
function jacobian_finitediff!(J, f, x, jac_config, integrator)
    (FiniteDiff.finite_difference_jacobian!(J, f, x, jac_config,
        dir = diffdir(integrator));
    2 * maximum(jac_config.colorvec))
end

function jacobian!(J::AbstractMatrix{<:Number}, f, x::AbstractArray{<:Number},
        fx::AbstractArray{<:Number}, integrator::DiffEqBase.DEIntegrator,
        jac_config)
    alg = unwrap_alg(integrator, true)
    if alg_autodiff(alg) isa AutoForwardDiff
        if integrator.iter == 1
            try
                forwarddiff_color_jacobian!(J, f, x, jac_config)
            catch e
                throw(FirstAutodiffJacError(e))
            end
        else
            forwarddiff_color_jacobian!(J, f, x, jac_config)
        end
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, maximum(jac_config.colorvec))
    elseif alg_autodiff(alg) isa AutoFiniteDiff
        isforward = alg_difftype(alg) === Val{:forward}
        if isforward
            forwardcache = get_tmp_cache(integrator, alg, unwrap_cache(integrator, true))[2]
            f(forwardcache, x)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            tmp = jacobian_finitediff_forward!(J, f, x, jac_config, forwardcache,
                integrator)
        else # not forward difference
            tmp = jacobian_finitediff!(J, f, x, jac_config, integrator)
        end
        integrator.stats.nf += tmp
    else
        error("$alg_autodiff not yet supported in jacobian! function")
    end
    nothing
end

function build_jac_config(alg, f::F1, uf::F2, du1, uprev, u, tmp, du2) where {F1, F2}
    haslinsolve = hasfield(typeof(alg), :linsolve)

    if !DiffEqBase.has_jac(f) && # No Jacobian if has analytical solution
       (!DiffEqBase.has_Wfact_t(f)) &&
       ((concrete_jac(alg) === nothing && (!haslinsolve || (haslinsolve && # No Jacobian if linsolve doesn't want it
           (alg.linsolve === nothing || LinearSolve.needs_concrete_A(alg.linsolve))))) ||
        (concrete_jac(alg) !== nothing && concrete_jac(alg))) # Jacobian if explicitly asked for
        jac_prototype = f.jac_prototype

        if jac_prototype isa SparseMatrixCSC
            if f.mass_matrix isa UniformScaling
                idxs = diagind(jac_prototype)
                @. @view(jac_prototype[idxs]) = 1
            else
                idxs = findall(!iszero, f.mass_matrix)
                @. @view(jac_prototype[idxs]) = @view(f.mass_matrix[idxs])
            end
        end

        sparsity, colorvec = sparsity_colorvec(f, u)
        if alg_autodiff(alg) isa AutoForwardDiff
            _chunksize = get_chunksize(alg) === Val(0) ? nothing : get_chunksize(alg) # SparseDiffEq uses different convection...
            T = if standardtag(alg)
                typeof(ForwardDiff.Tag(OrdinaryDiffEqTag(), eltype(u)))
            else
                typeof(ForwardDiff.Tag(uf, eltype(u)))
            end

            if _chunksize === Val{nothing}()
                _chunksize = nothing
            end
            jac_config = ForwardColorJacCache(uf, uprev, _chunksize; colorvec = colorvec,
                sparsity = sparsity, tag = T)
        elseif alg_autodiff(alg) isa AutoFiniteDiff
            if alg_difftype(alg) !== Val{:complex}
                jac_config = FiniteDiff.JacobianCache(tmp, du1, du2, alg_difftype(alg),
                    colorvec = colorvec,
                    sparsity = sparsity)
            else
                jac_config = FiniteDiff.JacobianCache(Complex{eltype(tmp)}.(tmp),
                    Complex{eltype(du1)}.(du1), nothing,
                    alg_difftype(alg), eltype(u),
                    colorvec = colorvec,
                    sparsity = sparsity)
            end
        else
            error("$alg_autodiff not yet supported in build_jac_config function")
        end
    else
        jac_config = nothing
    end
    jac_config
end

function get_chunksize(jac_config::ForwardDiff.JacobianConfig{
        T,
        V,
        N,
        D
}) where {T, V, N, D
}
    Val(N)
end # don't degrade compile time information to runtime information

function resize_jac_config!(jac_config::SparseDiffTools.ForwardColorJacCache, i)
    resize!(jac_config.fx, i)
    resize!(jac_config.dx, i)
    resize!(jac_config.t, i)
    ps = SparseDiffTools.adapt.(DiffEqBase.parameterless_type(jac_config.dx),
        SparseDiffTools.generate_chunked_partials(jac_config.dx,
            1:length(jac_config.dx),
            Val(ForwardDiff.npartials(jac_config.t[1]))))
    resize!(jac_config.p, length(ps))
    jac_config.p .= ps
end

function resize_jac_config!(jac_config::FiniteDiff.JacobianCache, i)
    resize!(jac_config, i)
    jac_config
end

function resize_grad_config!(grad_config::AbstractArray, i)
    resize!(grad_config, i)
    grad_config
end

function resize_grad_config!(grad_config::ForwardDiff.DerivativeConfig, i)
    resize!(grad_config.duals, i)
    grad_config
end

function resize_grad_config!(grad_config::FiniteDiff.GradientCache, i)
    @unpack fx, c1, c2 = grad_config
    fx !== nothing && resize!(fx, i)
    c1 !== nothing && resize!(c1, i)
    c2 !== nothing && resize!(c2, i)
    grad_config
end

function build_grad_config(alg, f::F1, tf::F2, du1, t) where {F1, F2}
    if !DiffEqBase.has_tgrad(f)
        if alg_autodiff(alg) isa AutoForwardDiff
            T = if standardtag(alg)
                typeof(ForwardDiff.Tag(OrdinaryDiffEqTag(), eltype(du1)))
            else
                typeof(ForwardDiff.Tag(f, eltype(du1)))
            end

            if du1 isa Array
                dualt = Dual{T, eltype(du1), 1}(first(du1) * t,
                    ForwardDiff.Partials((one(eltype(du1)),)))
                grad_config = similar(du1, typeof(dualt))
                fill!(grad_config, false)
            else
                grad_config = ArrayInterface.restructure(du1,
                    Dual{
                        T,
                        eltype(du1),
                        1
                    }.(du1,
                        (ForwardDiff.Partials((one(eltype(du1)),)),)) .*
                    false)
            end
        elseif alg_autodiff(alg) isa AutoFiniteDiff
            grad_config = FiniteDiff.GradientCache(du1, t, alg_difftype(alg))
        else
            error("$alg_autodiff not yet supported in build_grad_config function")
        end
    else
        grad_config = nothing
    end
    grad_config
end

function sparsity_colorvec(f, x)
    sparsity = f.sparsity

    if sparsity isa SparseMatrixCSC
        if f.mass_matrix isa UniformScaling
            idxs = diagind(sparsity)
            @. @view(sparsity[idxs]) = 1
        else
            idxs = findall(!iszero, f.mass_matrix)
            @. @view(sparsity[idxs]) = @view(f.mass_matrix[idxs])
        end
    end

    colorvec = DiffEqBase.has_colorvec(f) ? f.colorvec :
               (isnothing(sparsity) ? (1:length(x)) : matrix_colors(sparsity))
    sparsity, colorvec
end
