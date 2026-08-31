# Conventions follow qmat (https://github.com/Parallel-in-Time/qmat): nodes live
# on `[0, 1]` in increasing order, `Q[m, j] = ∫₀^{τₘ} ℓⱼ(s) ds` with `ℓⱼ` the
# Lagrange basis on the nodes, and `weights[j] = ∫₀¹ ℓⱼ(s) ds`.

"""
    SDCNodes

Distribution the collocation nodes are drawn from: `SDCNodes.Legendre` or
`SDCNodes.Equidistant`. Legendre nodes give the classical Gauss/Radau/Lobatto
collocation orders; equidistant nodes only reach the interpolatory order.
"""
@enumx SDCNodes begin
    Legendre
    Equidistant
end

"""
    SDCQuadrature

Which endpoints of the step are collocation nodes: `SDCQuadrature.Gauss`
(neither), `.RadauLeft`, `.RadauRight`, or `.Lobatto` (both). Only the two that
include the right endpoint support `SDCStepUpdate.LastNode`.
"""
@enumx SDCQuadrature begin
    Gauss
    RadauLeft
    RadauRight
    Lobatto
end

"""
    SDCSweeper

The sweep preconditioner `QΔ ≈ Q`, which sets how fast the iteration converges
but not what it converges to.

- `SDCSweeper.BE` — implicit Euler between the nodes, the original sweep.
- `SDCSweeper.FE` — explicit Euler between the nodes; strictly lower triangular,
  so the sweep needs no nonlinear solve at all.
- `SDCSweeper.Trapezoid` — the average of the two, second order.
- `SDCSweeper.LU` — Weiser's LU trick, `Uᵀ` from `Qᵀ = LU`. Makes the stiff-limit
  iteration matrix nilpotent and is usually the best serial choice.
- `SDCSweeper.Picard` — `QΔ = 0`, the unpreconditioned iteration.
- `SDCSweeper.BEpar` — `diag(τ)`, implicit Euler from the step start to each node.
- `SDCSweeper.MIN_SR_NS` — `diag(τ)/M`, which makes the non-stiff iteration
  matrix nilpotent, so the iteration terminates in `M` sweeps as `Δt → 0`.
- `SDCSweeper.MIN_SR_S` — tabulated entries that make the stiff-limit iteration
  matrix nilpotent instead, for stiff problems.
- `SDCSweeper.MIN_SR_FLEX` — `diag(τ)/k` on sweep `k`, changing every sweep, and
  falling back to `MIN_SR_S` past sweep `M`.

The last five are diagonal, so their sweeps decouple across the nodes. The three
`MIN_SR_*` families are from Čaklović, Lunet, Götschel and Ruprecht,
SIAM J. Sci. Comput. 47 (2025) A430-A453.
"""
@enumx SDCSweeper begin
    BE
    FE
    Trapezoid
    LU
    Picard
    BEpar
    MIN_SR_NS
    MIN_SR_S
    MIN_SR_FLEX
end

"""
    SDCStepUpdate

How the step solution is formed from the converged node values.
`SDCStepUpdate.Quadrature` applies the quadrature rule `uₙ + Δt Σ w_m f_m` and
works for any node set; `SDCStepUpdate.LastNode` takes `u_M`, which needs the
last node to be the right endpoint and gives up one order.
"""
@enumx SDCStepUpdate begin
    Quadrature
    LastNode
end

# Diagonal `QΔ`: the sweep decouples across the nodes.
const SDC_DIAGONAL_SWEEPERS = (
    SDCSweeper.Picard, SDCSweeper.BEpar, SDCSweeper.MIN_SR_NS,
    SDCSweeper.MIN_SR_S, SDCSweeper.MIN_SR_FLEX,
)

# Sweepers whose `QΔ` changes from one sweep to the next.
const SDC_SWEEP_DEPENDENT_SWEEPERS = (SDCSweeper.MIN_SR_FLEX,)

# Sweepers that are second order accurate on their own, and therefore gain two
# orders on the first sweep instead of one.
const SDC_SECOND_ORDER_SWEEPERS = (SDCSweeper.Trapezoid,)

# The monomial Vandermonde behind `Q` is ill-conditioned (cond ≈ 1e7 at M = 8),
# so coefficients are derived at this precision and rounded down at the end.
const SDC_COEFF_PRECISION = 256

struct SDCTableau{T}
    nodes::Vector{T}
    weights::Vector{T}
    Q::Matrix{T}
    # One entry per sweep for `MIN_SR_FLEX`, one entry in total otherwise.
    QΔ::Vector{Matrix{T}}
end

"""
    sdc_qdelta_for(tab, k)

The preconditioner for sweep `k`, counted from 1.
"""
@inline function sdc_qdelta_for(tab::SDCTableau, k::Int)
    return @inbounds tab.QΔ[min(k, length(tab.QΔ))]
end

"""
    _jacobi(n, α, β, x)

Evaluate the Jacobi polynomial ``P_n^{(α,β)}(x)`` by its three-term recurrence.
`α` and `β` are `0` or `1` here; `P_n^{(0,0)}` is the Legendre polynomial.
"""
function _jacobi(n::Int, α::Int, β::Int, x::T) where {T}
    n <= 0 && return one(T)
    p1 = (α + 1) + (α + β + 2) * (x - one(T)) / 2
    n == 1 && return p1
    p0 = one(T)
    for k in 2:n
        a = T(2k * (k + α + β) * (2k + α + β - 2))
        b = T((2k + α + β - 1) * (2k + α + β) * (2k + α + β - 2))
        c = T((2k + α + β - 1) * (α^2 - β^2))
        d = T(2 * (k + α - 1) * (k + β - 1) * (2k + α + β))
        p0, p1 = p1, ((b * x + c) * p1 - d * p0) / a
    end
    return p1
end

"""
    _bisect(f, a, b)

Bisect `f` on a bracket `[a, b]` with `f(a)f(b) < 0` down to the resolution of
the element type. Derivative free, so it is generic in the precision.

The iteration cap is not decoration. `BigFloat` has no smallest normal, so if a
bracket endpoint sits exactly on the root the midpoint halves towards it through
the whole exponent range and `m == a` never fires. Exact roots are removed by
the caller before bisection, and this bound catches anything that slips through.
"""
function _bisect(f::F, a::T, b::T) where {F, T}
    fa = f(a)
    iszero(fa) && return a
    positive = fa > 0
    for _ in 1:(8 * precision(T) + 64)
        m = (a + b) / 2
        (m == a || m == b) && return m
        fm = f(m)
        iszero(fm) && return m
        if (fm > 0) == positive
            a = m
        else
            b = m
        end
    end
    return (a + b) / 2
end

"""
    _jacobi_roots(n, α, β, T)

All `n` roots of ``P_n^{(α,β)}`` in `(-1, 1)`, in increasing order. The roots are
simple and, for the small `n` used by SDC, well separated, so a uniform
bracketing grid followed by bisection is both robust and generic in `T`.
"""
function _jacobi_roots(n::Int, α::Int, β::Int, ::Type{T}) where {T}
    n <= 0 && return T[]
    ngrid = 20n + 50
    roots = T[]
    xprev = -one(T)
    fprev = _jacobi(n, α, β, xprev)
    for i in 1:ngrid
        x = -one(T) + 2 * T(i) / T(ngrid)
        fx = _jacobi(n, α, β, x)
        if iszero(fx)
            # A grid point landed exactly on a root; this is the common case,
            # not a freak one — odd-degree Legendre polynomials have a root at
            # x = 0 and the grid is symmetric. Note it directly, and let the
            # `iszero(fprev)` guard below skip the following interval so it is
            # not counted twice. Bracketing it would hand `_bisect` an endpoint
            # that is exactly the root.
            push!(roots, x)
        elseif !iszero(fprev) && (fprev > 0) != (fx > 0)
            push!(roots, _bisect(y -> _jacobi(n, α, β, y), xprev, x))
        end
        xprev, fprev = x, fx
    end
    length(roots) == n || throw(
        ArgumentError(
            "SDC: found $(length(roots)) of $n roots for the Jacobi polynomial " *
                "P_$n^($α,$β); node generation failed"
        )
    )
    return roots
end

"""
    _sdc_nodes_big(M, node_type, quad_type)

The `M` collocation nodes on `[0, 1]`, in increasing order, at extended
precision. `quad_type` selects which interval endpoints are nodes:
`SDCQuadrature.Gauss` (neither), `.RadauLeft` (`0`), `.RadauRight` (`1`),
`.Lobatto` (both).
"""
function _sdc_nodes_big(M::Int, node_type::SDCNodes.T, quad_type::SDCQuadrature.T)
    T = BigFloat
    if node_type === SDCNodes.Equidistant
        # Closed form, matching qmat's EQUID convention.
        return if quad_type === SDCQuadrature.Gauss
            [T(m) / T(M + 1) for m in 1:M]
        elseif quad_type === SDCQuadrature.RadauLeft
            [T(m - 1) / T(M) for m in 1:M]
        elseif quad_type === SDCQuadrature.RadauRight
            [T(m) / T(M) for m in 1:M]
        else # SDCQuadrature.Lobatto
            [T(m - 1) / T(M - 1) for m in 1:M]
        end
    end
    # SDCNodes.Legendre — nodes on [-1, 1] from the Jacobi polynomials, then mapped.
    x = if quad_type === SDCQuadrature.Gauss
        _jacobi_roots(M, 0, 0, T)
    elseif quad_type === SDCQuadrature.RadauLeft
        vcat(-one(T), _jacobi_roots(M - 1, 0, 1, T))
    elseif quad_type === SDCQuadrature.RadauRight
        vcat(_jacobi_roots(M - 1, 1, 0, T), one(T))
    else # SDCQuadrature.Lobatto
        vcat(-one(T), _jacobi_roots(M - 2, 1, 1, T), one(T))
    end
    return (x .+ 1) ./ 2
end

"""
    _lagrange_integrals(τ, uppers)

`P[m, j] = ∫₀^{uppers[m]} ℓⱼ(s) ds` for the Lagrange basis `ℓⱼ` on the nodes `τ`.

Written in the monomial basis: the Lagrange coefficients are `V⁻¹` for the
Vandermonde matrix `V[k, i] = τₖ^{i-1}` (since `ℓⱼ(τₖ) = δⱼₖ`), and the monomial
integrals are `C[m, i] = uppers[m]^i / i`, so the result is `C V⁻¹`.
"""
function _lagrange_integrals(τ::Vector{BigFloat}, uppers::Vector{BigFloat})
    M = length(τ)
    V = [τ[k]^(i - 1) for k in 1:M, i in 1:M]
    C = [uppers[m]^i / i for m in 1:length(uppers), i in 1:M]
    return C / V
end

"""
    sdc_qdelta(T, sweeper, τ, Q)

The sweep preconditioner `QΔ ≈ Q`. Must be lower triangular for the sweep to
decouple into `M` successive `N`-sized solves; a diagonal `QΔ` decouples them
completely, which is what parallel-across-the-nodes SDC exploits.
"""
function sdc_qdelta(
        ::Type{T}, sweeper::SDCSweeper.T, τ::Vector{BigFloat}, Q::Matrix{BigFloat},
        node_type::SDCNodes.T, quad_type::SDCQuadrature.T, sweep::Int
    ) where {T}
    M = length(τ)
    QΔ = zeros(BigFloat, M, M)
    # Node spacings, with the first measured from the start of the step.
    δ = [i == 1 ? τ[1] : τ[i] - τ[i - 1] for i in 1:M]
    if sweeper === SDCSweeper.BE
        for i in 1:M, j in 1:i
            QΔ[i, j] = δ[j]
        end
    elseif sweeper === SDCSweeper.FE
        for i in 1:M, j in 1:(i - 1)
            QΔ[i, j] = δ[j + 1]
        end
    elseif sweeper === SDCSweeper.Trapezoid
        for i in 1:M
            for j in 1:i
                QΔ[i, j] += δ[j]
            end
            for j in 1:(i - 1)
                QΔ[i, j] += δ[j + 1]
            end
        end
        QΔ ./= 2
    elseif sweeper === SDCSweeper.LU
        # `check = false` because Q is singular whenever τ₁ = 0 (Lobatto and
        # Radau-left), where its first row vanishes. The resulting zero pivot is
        # the right answer: it makes the first node explicit, and u(τ₁) = uₙ.
        QΔ = Matrix(transpose(LinearAlgebra.lu(transpose(Q); check = false).U))
    elseif sweeper === SDCSweeper.Picard
        # QΔ stays zero, giving the unpreconditioned Picard iteration.
    elseif sweeper === SDCSweeper.BEpar
        for i in 1:M
            QΔ[i, i] = τ[i]
        end
    elseif sweeper === SDCSweeper.MIN_SR_NS
        for i in 1:M
            QΔ[i, i] = τ[i] / M
        end
    elseif sweeper === SDCSweeper.MIN_SR_S
        for (i, d) in enumerate(min_sr_s_coefficients(node_type, quad_type, M))
            QΔ[i, i] = d
        end
    elseif sweeper === SDCSweeper.MIN_SR_FLEX
        # Čaklović et al. Theorem 2.13: diag(τ)/k for the k-th sweep drives the
        # product of the stiff-limit iteration matrices to exactly zero after M
        # sweeps. Past that the theorem says nothing, so fall back to MIN-SR-S.
        if sweep <= M
            for i in 1:M
                QΔ[i, i] = τ[i] / sweep
            end
        else
            for (i, d) in enumerate(min_sr_s_coefficients(node_type, quad_type, M))
                QΔ[i, i] = d
            end
        end
    else
        throw(ArgumentError("SDC: unknown sweeper $(sweeper)"))
    end
    return T.(QΔ)
end

"""
    SDCTableau(T, M, node_type, quad_type, sweeper)

Build every coefficient array the sweep needs, in the element type `T`.
"""
function SDCTableau(
        ::Type{T}, M::Int, node_type::SDCNodes.T, quad_type::SDCQuadrature.T,
        sweeper::SDCSweeper.T, num_sweeps::Int = 1
    ) where {T}
    return setprecision(BigFloat, SDC_COEFF_PRECISION) do
        τ = _sdc_nodes_big(M, node_type, quad_type)
        Q = _lagrange_integrals(τ, τ)
        weights = vec(_lagrange_integrals(τ, [one(BigFloat)]))
        nsweeps = sweeper in SDC_SWEEP_DEPENDENT_SWEEPERS ? max(1, num_sweeps) : 1
        QΔ = [
            sdc_qdelta(T, sweeper, τ, Q, node_type, quad_type, k) for k in 1:nsweeps
        ]
        SDCTableau{T}(T.(τ), T.(weights), T.(Q), QΔ)
    end
end

"""
    min_sr_s_coefficients(node_type, quad_type, M)

Tabulated MIN-SR-S diagonal entries, or an error naming what is missing.
"""
function min_sr_s_coefficients(node_type::SDCNodes.T, quad_type::SDCQuadrature.T, M::Int)
    coeffs = get(MIN_SR_S_COEFFICIENTS, (node_type, quad_type, M), nothing)
    coeffs === nothing && throw(
        ArgumentError(
            "SDC: no tabulated MIN-SR-S coefficients for $(M) $(node_type)/$(quad_type) " *
                "nodes; the available node counts are $(sort(unique(k[3] for k in keys(MIN_SR_S_COEFFICIENTS))))"
        )
    )
    return coeffs
end

"""
    sdc_collocation_order(M, node_type, quad_type)

Order of the underlying collocation method, which caps the order SDC can reach
no matter how many sweeps are taken.

Legendre nodes give the classical Gauss/Radau/Lobatto orders `2M`, `2M-1`,
`2M-2`. Equidistant nodes give the interpolatory order `M`, raised to `M+1` for
even `M` on the symmetric (Gauss and Lobatto) rules. Matches
`qmat.qcoeff.collocation.Collocation.order`.
"""
function sdc_collocation_order(M::Int, node_type::SDCNodes.T, quad_type::SDCQuadrature.T)
    if node_type === SDCNodes.Legendre
        return if quad_type === SDCQuadrature.Gauss
            2M
        elseif quad_type === SDCQuadrature.Lobatto
            2M - 2
        else # RadauLeft, RadauRight
            2M - 1
        end
    end
    if quad_type === SDCQuadrature.Gauss || quad_type === SDCQuadrature.Lobatto
        return M + (M % 2)
    end
    return M
end
