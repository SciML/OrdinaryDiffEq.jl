using OrdinaryDiffEqFIRK, LinearSolve, LinearAlgebra, Test
using SciMLOperators: AbstractSciMLOperator

const N_brus = 8
const xyd_brus = range(0, stop = 1, length = N_brus)

brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
limit_brus(a, N) = a == N + 1 ? 1 : a == 0 ? N : a

function brusselator_2d_loop!(du, u, p, t)
    A, B, alpha, dx = p
    alpha = alpha / dx^2
    @inbounds for I in CartesianIndices((N_brus, N_brus))
        i, j = Tuple(I)
        x, y = xyd_brus[i], xyd_brus[j]
        ip1 = limit_brus(i + 1, N_brus)
        im1 = limit_brus(i - 1, N_brus)
        jp1 = limit_brus(j + 1, N_brus)
        jm1 = limit_brus(j - 1, N_brus)
        du[i, j, 1] = alpha *
            (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] - 4u[i, j, 1]) +
            B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] +
            brusselator_f(x, y, t)
        du[i, j, 2] = alpha *
            (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] - 4u[i, j, 2]) +
            A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
    end
    return nothing
end

function init_brusselator_2d(xyd)
    N = length(xyd)
    u = zeros(N, N, 2)
    for I in CartesianIndices((N, N))
        x = xyd[I[1]]
        y = xyd[I[2]]
        u[I, 1] = 22 * (y * (1 - y))^(3 / 2)
        u[I, 2] = 27 * (x * (1 - x))^(3 / 2)
    end
    return u
end

brus_prob = ODEProblem(
    brusselator_2d_loop!, init_brusselator_2d(xyd_brus), (0.0, 1.0),
    (3.4, 1.0, 10.0, step(xyd_brus))
)

@testset "RadauIIA5 matrix-free Krylov" begin
    integ = init(brus_prob, RadauIIA5(linsolve = KrylovJL_GMRES()))
    @test integ.cache.J isa AbstractSciMLOperator
    @test integ.cache.W2 isa OrdinaryDiffEqFIRK.ComplexShiftOperator

    ref = solve(brus_prob, RadauIIA5(), abstol = 1.0e-12, reltol = 1.0e-12)
    sol = solve(
        brus_prob, RadauIIA5(linsolve = KrylovJL_GMRES()),
        abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ ref.u[end] rtol = 1.0e-7

    sol_bicg = solve(
        brus_prob, RadauIIA5(linsolve = KrylovJL_BICGSTAB()),
        abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(sol_bicg)
    @test sol_bicg.u[end] ≈ ref.u[end] rtol = 1.0e-7
end

@testset "RadauIIA5 matrix-free Krylov with mass matrix" begin
    n_heat = 24
    function heat_nonlinear!(du, u, p, t)
        dx2 = inv((1 / (length(u) + 1))^2)
        @inbounds for i in eachindex(u)
            um = i == 1 ? zero(eltype(u)) : u[i - 1]
            up = i == length(u) ? zero(eltype(u)) : u[i + 1]
            du[i] = dx2 * (um - 2u[i] + up) - u[i]^2
        end
        return nothing
    end
    u0 = [sinpi(i / (n_heat + 1)) for i in 1:n_heat]
    mm = Diagonal(range(1.0, 2.0, length = n_heat))
    prob = ODEProblem(
        ODEFunction(heat_nonlinear!, mass_matrix = mm), u0, (0.0, 1.0)
    )

    ref = solve(prob, RadauIIA5(), abstol = 1.0e-12, reltol = 1.0e-12)
    sol = solve(
        prob, RadauIIA5(linsolve = KrylovJL_GMRES()), abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ ref.u[end] rtol = 1.0e-7
end

@testset "Unsupported FIRK methods report what is supported" begin
    for alg in (
            RadauIIA3(linsolve = KrylovJL_GMRES()),
            RadauIIA9(linsolve = KrylovJL_GMRES()),
            AdaptiveRadau(linsolve = KrylovJL_GMRES()),
            GaussLegendre(linsolve = KrylovJL_GMRES()),
        )
        @test_throws "RadauIIA5" init(brus_prob, alg)
    end
end
