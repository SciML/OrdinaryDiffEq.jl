using Test
using DiffEqBase
using DynamicQuantities
using LinearAlgebra

@testset "DiffEqBaseDynamicQuantitiesExt" begin
    # Basic quantity hooks
    q = 3.0u"m"
    @test DiffEqBase.ODE_DEFAULT_NORM(q, 0.0) == 3.0

    qc = (3.0 + 4.0im)u"m"
    @test DiffEqBase.UNITLESS_ABS2(qc) == 25.0

    r = DiffEqBase._rate_prototype(2.0u"m", 4.0u"s", 1)
    @test isapprox(ustrip(r), 2.0)
    @test oneunit(r) == oneunit(1.0u"m") / oneunit(1.0u"s")

    dt = DiffEqBase.timedepentdtmin(1.0u"s", 1.0u"ms")
    @test isapprox(ustrip(dt), 0.001)
    @test oneunit(dt) == oneunit(1.0u"s")

    # Factorization bridge for Quantity matrices
    A = [2.0u"m" 0.0u"m"; 0.0u"m" 4.0u"m"]
    b = [4.0u"m", 8.0u"m"]

    W = DiffEqBase.default_factorize(A)
    x = W \ b
    @test maximum(abs.(ustrip.(A * x .- b))) ≤ 1.0e-12

    x2 = similar(b)
    ldiv!(x2, W, b)
    @test maximum(abs.(ustrip.(A * x2 .- b))) ≤ 1.0e-12

    # _infer_ut fallback when all entries are zero
    Az = fill(0.0u"m", 2, 2)
    Wz = DiffEqBase.default_factorize(Az)
    @test Wz.ut == oneunit(1.0)

    # empty-matrix path (coverage + sanity)
    A0 = Matrix{typeof(1.0u"m")}(undef, 0, 0)
    W0 = DiffEqBase.default_factorize(A0)
    b0 = typeof(1.0u"m")[]
    @test length(W0 \ b0) == 0
end
