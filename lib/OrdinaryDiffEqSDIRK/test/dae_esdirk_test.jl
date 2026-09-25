using OrdinaryDiffEqSDIRK, LinearAlgebra, Test
using SciMLBase: alg_order
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit, ShampineCollocationInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff
import DifferentiationInterface as DI

afd_cs3 = AutoForwardDiff(chunksize = 3)
function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end
function rober(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    return [
        -k₁ * y₁ + k₃ * y₂ * y₃,
        k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
        y₁ + y₂ + y₃ - 1,
    ]
end
M = Diagonal([1.0, 1.0, 0.0])
roberf = ODEFunction{true, SciMLBase.AutoSpecialize}(rober, mass_matrix = M)
roberf_oop = ODEFunction{false, SciMLBase.AutoSpecialize}(rober, mass_matrix = M)
prob_mm = ODEProblem(roberf, [1.0, 0.0, 0.2], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
prob_mm_oop = ODEProblem(roberf_oop, [1.0, 0.0, 0.2], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
# Both should be inferable so long as AutoSpecialize is used...
sol = @inferred solve(
    prob_mm, KenCarp47(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)
sol = @inferred solve(
    prob_mm_oop, KenCarp47(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)

# These tests flex differentiation of the solver and through the initialization
# To only test the solver part and isolate potential issues, set the initialization to consistent
@testset "Inplace: $(isinplace(_prob)), BrownBasic: $(initalg isa BrownFullBasicInit), Autodiff: $autodiff" for _prob in [
            prob_mm, prob_mm_oop,
        ],
        initalg in [BrownFullBasicInit(), ShampineCollocationInit()],
        autodiff in [AutoForwardDiff(chunksize = 3), AutoFiniteDiff()]

    alg = KenCarp47(; autodiff)
    function f(p)
        sol = @inferred solve(
            remake(_prob; p), alg, abstol = 1.0e-14,
            reltol = 1.0e-14, initializealg = initalg
        )
        sum(sol)
    end
    @test DI.gradient(f, AutoForwardDiff(), [0.04, 3.0e7, 1.0e4]) ≈ [0, 0, 0] atol = 1.0e-8
end

# The rober mass matrix above is `Diagonal([1, 1, 0])`, where plain `f` and `M⁻¹f`
# happen to agree: the differential rows have `Mᵢᵢ = 1`, and the algebraic row has
# `f = 0` at a consistent point. That hides a wrong FSAL <-> stage-variable
# conversion, which costs an `O(dt)` blunder on the first step and pins the global
# order at 1. A non-unit diagonal `M` separates the two conventions on every row.
const SDIRK = OrdinaryDiffEqSDIRK

# `f = M * (A * u)` makes `M u' = f` exactly equivalent to `u' = A * u`, so the
# mass-matrix handling is the only thing under test and the reference is exact.
const A_lin = [-0.5 0.1; 0.2 -0.3]
const u0_lin = [1.0, -0.5]
const T_lin = 1.0
const exact_lin = exp(A_lin * T_lin) * u0_lin

function lin_mm_prob(M; inplace = true)
    f = if inplace
        (du, u, p, t) -> (mul!(du, A_lin, u); du .= M * du; nothing)
    else
        (u, p, t) -> M * (A_lin * u)
    end
    return ODEProblem(ODEFunction{inplace}(f; mass_matrix = M), u0_lin, (0.0, T_lin))
end

function observed_order(prob, alg)
    errs = [
        norm(solve(prob, alg; dt = 1 / 2^k, adaptive = false).u[end] - exact_lin)
            for k in 5:8
    ]
    return maximum(log2(errs[i] / errs[i + 1]) for i in 1:(length(errs) - 1))
end

# Explicit-first-stage methods recover z₁ from the FSAL vector — the path that
# needs the conversion.
esdirk_algs = [Kvaerno3(), KenCarp3(), Kvaerno4(), KenCarp4(), TRBDF2(), ESDIRK325L2SA()]

@testset "Non-unit diagonal mass matrix keeps full order" begin
    M = Diagonal([2.0, 3.0])
    @testset "$(nameof(typeof(alg))), inplace=$ip" for alg in esdirk_algs, ip in (true, false)
        @test observed_order(lin_mm_prob(M; inplace = ip), alg) ≥ alg_order(alg) - 0.35
    end
end

@testset "Identity mass matrix is unaffected" begin
    @testset "$(nameof(typeof(alg)))" for alg in esdirk_algs
        @test observed_order(lin_mm_prob(Diagonal([1.0, 1.0])), alg) ≥ alg_order(alg) - 0.35
    end
end

@testset "Non-diagonal mass matrix" begin
    prob = lin_mm_prob([2.0 1.0; 1.0 3.0])
    # Restricted: the elementwise M⁻¹ in the explicit first stage cannot represent it.
    for alg in esdirk_algs
        @test_throws ArgumentError solve(prob, alg; dt = 0.01, adaptive = false)
    end
    # Implicit-first-stage SDIRKs solve for z₁ and keep general mass-matrix support.
    for alg in (Cash4(), Hairer4(), SDIRK2())
        @test !SDIRK.only_diagonal_mass_matrix(alg)
        @test observed_order(prob, alg) ≥ alg_order(alg) - 0.35
    end
end

@testset "only_diagonal_mass_matrix matches explicit_first_stage" begin
    # The trait is hand-listed so it stays constant-foldable in the hot path; this
    # keeps it from drifting away from the tableau flag that motivates it.
    all_algs = [
        ImplicitEuler(), ImplicitMidpoint(), Trapezoid(), TRBDF2(), SDIRK2(),
        Kvaerno3(), KenCarp3(), Cash4(), Hairer4(), Hairer42(), SSPSDIRK2(),
        Kvaerno4(), Kvaerno5(), KenCarp4(), KenCarp47(), KenCarp5(), KenCarp58(),
        ESDIRK54I8L2SA(), SFSDIRK4(), SFSDIRK5(), SFSDIRK6(), SFSDIRK7(), SFSDIRK8(),
        ESDIRK325L2SA(), ESDIRK436L2SA2(), ESDIRK437L2SA(), ESDIRK547L2SA2(),
        ESDIRK659L2SA(), CFNLIRK3(), ARS343(), ARS222(), ARS232(), ARS443(),
        IMEXSSP222(), IMEXSSP2322(), IMEXSSP3332(), IMEXSSP3433(), BHR553(),
    ]
    for alg in all_algs
        tab = SDIRK.ESDIRKIMEXTableau(alg, Float64, Float64)
        @testset "$(nameof(typeof(alg)))" begin
            @test SDIRK.only_diagonal_mass_matrix(alg) == tab.explicit_first_stage
        end
    end
end
