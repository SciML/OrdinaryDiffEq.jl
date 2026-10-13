using OrdinaryDiffEq, OrdinaryDiffEqSDIRK, Test

decay!(du, u, p, t) = (du .= -p .* u; nothing)

const DELETEAT_ADDAT_ALGS = (
    Trapezoid(), ImplicitEuler(), KenCarp4(), Tsit5(), Vern7(), Rodas5P(),
)

function reduced_reference(integ, alg, p, tend; kwargs...)
    prob = ODEProblem(decay!, copy(integ.u), (integ.t, tend), p)
    return solve(prob, alg; kwargs...).u[end]
end

@testset "deleteat!/addat! mid-solve with $(nameof(typeof(alg)))" for alg in
    DELETEAT_ADDAT_ALGS
    tend = 1.0
    kwargs = (; adaptive = false, dt = 1 / 64, save_everystep = false)

    @testset "deleteat!" begin
        p = collect(1.0:6.0)
        integ = init(ODEProblem(decay!, collect(1.0:6.0), (0.0, tend), p), alg; kwargs...)
        for _ in 1:5
            step!(integ)
        end
        kept = integ.u[[1, 4, 5, 6]]
        deleteat!(integ, [2, 3])
        deleteat!(p, [2, 3])
        @test integ.u == kept
        ref = reduced_reference(integ, alg, copy(p), tend; kwargs...)
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test integ.t == tend
        @test integ.u ≈ ref rtol = 1.0e-12
    end

    @testset "addat!" begin
        p = collect(1.0:6.0)
        integ = init(ODEProblem(decay!, collect(1.0:6.0), (0.0, tend), p), alg; kwargs...)
        for _ in 1:5
            step!(integ)
        end
        old = copy(integ.u)
        addat!(integ, 7:8)
        append!(p, [7.0, 8.0])
        @test integ.u[1:6] == old
        integ.u[7:8] .= [10.0, 20.0]
        u_modified!(integ, true)
        ref = reduced_reference(integ, alg, copy(p), tend; kwargs...)
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test integ.t == tend
        @test integ.u ≈ ref rtol = 1.0e-12
    end
end
