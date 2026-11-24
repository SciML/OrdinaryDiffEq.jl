# Evaluate F_k = [f]_k and J_k = [f_u]_k along the jet
# u(σ) = u + Σ_{j=1}^{Q-1} v_j σ^j, t(σ) = t + hσ. The v_j are free stage
# inputs, with ∂F_k/∂v_j = J_{k-j}. The last column of V is outside the jet.
# The implicit solver requests only J_0 for its simplified Newton matrix.
struct StageBundle{F}
    f::F
    constant_jacobian::Bool   # every requested J_k is independent of (u, V, t, h): exact reuse
end
@inline function (bundle::StageBundle)(F, u, V, t, h)
    bundle.f(PolynomialOutputs((F,)), u, V, t, h)
    return nothing
end
@inline function (bundle::StageBundle)(F, J, u, V, t, h)
    bundle.f(PolynomialOutputs((F, J)), u, V, t, h)
    return nothing
end

function build_stage_bundles(
        f::ODEFunction{iip}, p, Q::Int, length::Int;
        jacobian_order::Int = Q - 1
    ) where {iip}
    return build_stage_bundles(unwrapped_f(f), Val(iip), p, Q, length; jacobian_order)
end

# In-place kernels values(F, u, V, t, h) and jacobians(F, J, u, V, t, h).
# F is n×Q and J is n×n×(jacobian_order + 1).
function build_stage_bundles(
        f, ::Val{iip}, p, Q::Int, length::Int;
        jacobian_order::Int = Q - 1
    ) where {iip}
    Q >= 1 || throw(ArgumentError("stage kernels need Q ≥ 1"))
    0 <= jacobian_order < Q || throw(ArgumentError("jacobian_order must satisfy 0 ≤ jacobian_order < Q"))
    kernels = get!(
        Dict{Tuple, Vector{PolynomialBundleCacheEntry}},
        POLYNOMIAL_BUNDLE_CACHE, f
    )::Dict{Tuple, Vector{PolynomialBundleCacheEntry}}
    key = (:stages, Q, length, iip, jacobian_order)
    entries = get!(Vector{PolynomialBundleCacheEntry}, kernels, key)
    index = findfirst(entry -> entry.parameters == p, entries)
    index === nothing || return entries[index].value
    @variables t0::Real dt::Real
    u0 = Symbolics.variables(:u0, 1:length)
    V = Symbolics.variables(:v, 1:length, 1:Q)
    if Q == 1
        # Order-zero jet: plain evaluation at (u0, t0).
        if iip
            f0 = similar(u0)
            f(f0, u0, p, t0)
        else
            f0 = f(u0, p, t0)
        end
        coefficients = [[f0[i]] for i in 1:length]
    else
        jet = [TaylorScalar(ntuple(k -> k == 1 ? u0[i] : V[i, k - 1], Q)) for i in 1:length]
        tjet = TaylorScalar{Q - 1}(t0, dt)
        if iip
            fu = similar(jet)
            f(fu, jet, p, tjet)
        else
            fu = f(jet, p, tjet)
        end
        coefficients = [collect(TaylorDiff.flatten(fu[i])) for i in 1:length]
    end
    Fmat = [coefficients[i][k + 1] for i in 1:length, k in 0:(Q - 1)]
    Jstack = cat([Symbolics.jacobian(Fmat[:, k + 1], u0) for k in 0:jacobian_order]...; dims = 3)
    _, values! = build_function(
        [Fmat], u0, V, t0, dt;
        expression = Val(false), cse = true, skipzeros = false
    )
    _, jacobians! = build_function(
        [Fmat, Jstack], u0, V, t0, dt;
        expression = Val(false), cse = true, skipzeros = false
    )
    constant_jacobian = all(x -> isempty(Symbolics.get_variables(x)), Jstack)
    bundles = (StageBundle(values!, constant_jacobian), StageBundle(jacobians!, constant_jacobian))
    push!(entries, PolynomialBundleCacheEntry(p, bundles))
    return bundles
end
