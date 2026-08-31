using OrdinaryDiffEqSDC
using OrdinaryDiffEqSDC: SDCTableau, sdc_collocation_order, sdc_solver_index
using LinearAlgebra
using Test

include("sdc_reference.jl")

@testset "SDC coefficients match qmat" begin
    for case in QMAT_GOLDEN
        label = "M=$(case.M) $(case.node_type)/$(case.quad_type)"
        for (sweeper, QDelta) in case.QDelta
            tab = SDCTableau(
                Float64, case.M, case.node_type, case.quad_type, sweeper
            )
            @testset "$label $sweeper" begin
                @test tab.nodes ≈ case.nodes atol = 1.0e-13
                @test tab.weights ≈ case.weights atol = 1.0e-13
                @test tab.Q ≈ case.Q atol = 1.0e-13
                @test first(tab.QΔ) ≈ QDelta atol = 1.0e-13
            end
        end
        @test sdc_collocation_order(case.M, case.node_type, case.quad_type) ==
            case.coll_order
    end
end

@testset "SDC coefficient invariants" begin
    for node_type in instances(SDCNodes.T), quad_type in instances(SDCQuadrature.T)
        minimum_nodes = quad_type in (SDCQuadrature.Lobatto, SDCQuadrature.RadauLeft) ? 2 : 1
        for M in max(2, minimum_nodes):6
            tab = SDCTableau(Float64, M, node_type, quad_type, SDCSweeper.BE)
            τ = tab.nodes
            @testset "$node_type/$quad_type M=$M" begin
                @test issorted(τ)
                @test all(0 .<= τ .<= 1)
                quad_type in (SDCQuadrature.Lobatto, SDCQuadrature.RadauLeft) && @test τ[1] == 0
                quad_type in (SDCQuadrature.Lobatto, SDCQuadrature.RadauRight) && @test τ[end] == 1
                # Q ⋅ τ^{n} integrates monomials exactly up to degree M-1.
                for n in 0:(M - 1)
                    @test tab.Q * (τ .^ n) ≈ (τ .^ (n + 1)) ./ (n + 1) atol = 1.0e-13
                    @test sum(tab.weights .* τ .^ n) ≈ 1 / (n + 1) atol = 1.0e-13
                end
                @test first(tab.QΔ) ≈ LowerTriangular(first(tab.QΔ))
                @test first(tab.QΔ) * ones(M) ≈ τ atol = 1.0e-14
            end
        end
    end
end

@testset "SDC sweeper shapes" begin
    M = 4
    for sweeper in instances(SDCSweeper.T)
        tab = SDCTableau(Float64, M, SDCNodes.Legendre, SDCQuadrature.RadauRight, sweeper)
        @test first(tab.QΔ) ≈ LowerTriangular(first(tab.QΔ))
        # The sweepers advertised as parallel-ready must actually be diagonal.
        if sweeper in OrdinaryDiffEqSDC.SDC_DIAGONAL_SWEEPERS
            @test first(tab.QΔ) ≈ Diagonal(first(tab.QΔ))
        else
            @test !(first(tab.QΔ) ≈ Diagonal(first(tab.QΔ)))
        end
        if sweeper === SDCSweeper.FE
            @test all(iszero, diag(first(tab.QΔ)))
        end
    end
    be = SDCTableau(Float64, M, SDCNodes.Legendre, SDCQuadrature.RadauRight, SDCSweeper.BE).QΔ |> first
    fe = SDCTableau(Float64, M, SDCNodes.Legendre, SDCQuadrature.RadauRight, SDCSweeper.FE).QΔ |> first
    tr = SDCTableau(Float64, M, SDCNodes.Legendre, SDCQuadrature.RadauRight, SDCSweeper.Trapezoid).QΔ |> first
    @test tr ≈ (be .+ fe) ./ 2
end

@testset "SDC explicit-node bookkeeping" begin
    # τ₁ = 0 for Lobatto and Radau-left, so node 1 has no implicit solve.
    for quad_type in (SDCQuadrature.Lobatto, SDCQuadrature.RadauLeft), sweeper in (SDCSweeper.BE, SDCSweeper.BEpar, SDCSweeper.MIN_SR_NS)
        tab = SDCTableau(Float64, 4, SDCNodes.Legendre, quad_type, sweeper)
        index = sdc_solver_index(tab.QΔ)
        @test index[1] == 0
        @test index[2:end] == collect(1:3)
    end
    tab = SDCTableau(Float64, 4, SDCNodes.Legendre, SDCQuadrature.RadauRight, SDCSweeper.BE)
    @test sdc_solver_index(tab.QΔ) == collect(1:4)
    tab = SDCTableau(Float64, 4, SDCNodes.Legendre, SDCQuadrature.RadauRight, SDCSweeper.FE)
    @test all(iszero, sdc_solver_index(tab.QΔ))
end

@testset "SDC coefficients are precision generic" begin
    small = SDCTableau(Float64, 5, SDCNodes.Legendre, SDCQuadrature.RadauRight, SDCSweeper.LU)
    big = SDCTableau(BigFloat, 5, SDCNodes.Legendre, SDCQuadrature.RadauRight, SDCSweeper.LU)
    @test Float64.(big.nodes) ≈ small.nodes atol = 1.0e-15
    @test Float64.(big.Q) ≈ small.Q atol = 1.0e-15
    τ = big.nodes
    # Extended precision exposes any conditioning loss in the Q construction.
    @test maximum(abs, big.Q * (τ .^ 3) - (τ .^ 4) ./ 4) < 1.0e-40
end

@testset "SDC argument validation" begin
    # Unknown names are rejected by the enum types themselves.
    @test_throws Exception SDC(node_type = :Chebyshev)
    @test_throws Exception SDC(sweeper = :NotASweeper)
    @test_throws ArgumentError SDC(num_sweeps = -1)
    @test_throws ArgumentError SDC(num_nodes = 1, quad_type = SDCQuadrature.Lobatto)
    @test_throws ArgumentError SDC(
        num_sweeps = 0, step_update = SDCStepUpdate.LastNode
    )
    @test_throws ArgumentError SDC(quad_type = SDCQuadrature.Gauss, step_update = SDCStepUpdate.LastNode)
    @test SDC(quad_type = SDCQuadrature.RadauRight, step_update = SDCStepUpdate.LastNode) isa SDC
end

@testset "SDC alg_order" begin
    for K in 0:8
        alg = SDC(num_nodes = 3, quad_type = SDCQuadrature.RadauRight, num_sweeps = K)
        @test OrdinaryDiffEqSDC.alg_order(alg) == max(1, min(K + 1, 5))
    end
    @test OrdinaryDiffEqSDC.alg_order(
        SDC(num_nodes = 4, quad_type = SDCQuadrature.Gauss, num_sweeps = 1, sweeper = SDCSweeper.Trapezoid)
    ) == 3
    @test OrdinaryDiffEqSDC.alg_order(
        SDC(num_nodes = 4, quad_type = SDCQuadrature.Gauss, num_sweeps = 1, sweeper = SDCSweeper.BE)
    ) == 2
    @test OrdinaryDiffEqSDC.alg_order(
        SDC(num_nodes = 4, quad_type = SDCQuadrature.RadauRight, num_sweeps = 3, step_update = SDCStepUpdate.LastNode)
    ) == 3
end
