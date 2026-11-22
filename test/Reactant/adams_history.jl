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
    function endpoint(u, p, nsteps)
        prob = ODEProblem(rhs, u, (32.0, 32.0 + 3.5 * nsteps), p)
        integ = init(prob, alg; dt = 3.5, dtmax = 3.5, adaptive = false, save_everystep = false)
        sol = solve!(integ)
        lags = order == 3 ? (integ.cache.k2, integ.cache.k3) : (integ.cache.k2, integ.cache.k3, integ.cache.k4)
        return (; sol.retcode, integ.iter, integ.cache.step, u = sol.u[end], lags, integ.fsalfirst)
    end

    u, p = [1 / 64], [Float64(rate)]
    host = endpoint(u, p, 5)
    compiled = Reactant.@jit endpoint(Reactant.to_rarray(u), Reactant.to_rarray(p), 5)
    @test Bool(compiled.retcode == host.retcode == ReturnCode.Success)
    @test Int(compiled.iter) == host.iter == 5
    @test Int(compiled.step) == host.step == order
    @test only(host.u) ≈ Float64(trajectory[end]) rtol = 16eps()
    @test only(Array(compiled.u)) ≈ Float64(trajectory[end]) rtol = 16eps()
    for j in eachindex(host.lags)
        expected = Float64(rate * trajectory[end - j])
        @test only(host.lags[j]) ≈ expected rtol = 16eps()
        @test only(Array(compiled.lags[j])) ≈ expected rtol = 16eps()
    end
    @test only(host.fsalfirst) ≈ Float64(rate * trajectory[end - 1]) rtol = 16eps()
    @test only(Array(compiled.fsalfirst)) ≈ Float64(rate * trajectory[end - 1]) rtol = 16eps()

    if iip
        for nsteps in 1:(order - 1)
            host = endpoint(u, p, nsteps)
            compiled = Reactant.@jit endpoint(Reactant.to_rarray(u), Reactant.to_rarray(p), nsteps)
            for lag in 2:order
                j = order - lag + 1
                expected = j <= nsteps ? Float64(rate * trajectory[j]) : 0.0
                @test only(host.lags[lag - 1]) ≈ expected rtol = 16eps()
                @test only(Array(compiled.lags[lag - 1])) ≈ expected rtol = 16eps()
            end
        end
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
