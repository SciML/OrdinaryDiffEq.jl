using OrdinaryDiffEqAdamsBashforthMoulton: AB3, AB4, AB5, ABM32, ABM43
using SciMLBase: ODEProblem, init, solve!, ReturnCode
using Reactant
using Test

@testset "Compiled Adams history" for (alg, order) in ((AB3(), 3), (AB4(), 4)), iip in (false, true)
    rate = -1 // 8
    dt = 7 // 2
    z = rate * dt
    startup_order = order == 3 ? 2 : 4
    startup = sum(z^j / factorial(j) for j in 0:startup_order)
    weights = order == 3 ? (23, -16, 5) .// 12 : (55, -59, 37, -9) .// 24
    trajectory = [big(1) // 64 * startup^j for j in 0:(order - 1)]
    while length(trajectory) <= 5
        increment = sum(w * trajectory[end - j + 1] for (j, w) in enumerate(weights))
        push!(trajectory, last(trajectory) + z * increment)
    end

    rhs = iip ? ((du, u, p, t) -> (du .= p .* u; nothing)) : ((u, p, t) -> p .* u)
    function endpoint(u, p)
        prob = ODEProblem(rhs, u, (32.0, 49.5), p)
        integ = init(prob, alg; dt = 3.5, dtmax = 3.5, adaptive = false, save_everystep = false)
        sol = solve!(integ)
        return sol.retcode, integ.iter, integ.cache.step, sol.u[end], integ.cache.k2, integ.cache.k3
    end

    u, p = [1 / 64], [Float64(rate)]
    host = endpoint(u, p)
    compiled = Reactant.@jit endpoint(Reactant.to_rarray(u), Reactant.to_rarray(p))
    @test Bool(compiled[1] == host[1] == ReturnCode.Success)
    @test Int(compiled[2]) == host[2] == 5
    @test Int(compiled[3]) == host[3] == order
    for (j, expected) in enumerate((trajectory[end], rate * trajectory[end - 1], rate * trajectory[end - 2]))
        @test only(host[3 + j]) ≈ Float64(expected) rtol = 16eps()
        @test only(Array(compiled[3 + j])) ≈ Float64(expected) rtol = 16eps()
    end
end

@testset "Compiled Adams startup rejection" for alg in (ABM32(), ABM43(), AB5()), iip in (false, true)
    rhs = iip ? ((du, u, p, t) -> (du .= p .* u; nothing)) : ((u, p, t) -> p .* u)
    function setup(u, p)
        return init(ODEProblem(rhs, u, (32.0, 49.5), p), alg; dt = 3.5, adaptive = false)
    end
    u, p = [1 / 64], [-1 / 8]
    integ = setup(u, p)
    @test solve!(integ).retcode == ReturnCode.Success
    @test integ.iter == 5
    @test integ.cache.step == (alg isa ABM32 ? 3 : alg isa ABM43 ? 4 : 5)
    startup = (u, p) -> setup(u, p).cache.step
    @test_throws ArgumentError Reactant.@jit startup(Reactant.to_rarray(u), Reactant.to_rarray(p))
end
