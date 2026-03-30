# Fetch packages.
using DiffEqDevTools, NonlinearSolve, Plots, Test

let
    # Prepares NonlinearProblem.
    f(u, p) = 3u .^ 3 .+ 2u .^ 2 .+ u + .-p
    u0 = [1.0, 6.0]
    p = [1.0, 3.0]
    static_prob = NonlinearProblem(f, u0, p)
    real_sol = solve(static_prob, NewtonRaphson(), reltol = 1e-15, abstol = 1e-15)

    # Sets WP input.
    abstols = 1.0 ./ 10.0 .^ (8:12)
    reltols = 1.0 ./ 10.0 .^ (8:12)
    setups = [Dict(:alg => NewtonRaphson())
              Dict(:alg => TrustRegion())]
    solnames = ["NewtonRaphson"; "TrustRegion"]

    # Makes WP-diagram
    wp = WorkPrecisionSet(static_prob, abstols, reltols, setups; names = solnames,
        numruns = 100, appxsol = real_sol, error_estimate = :l2)

    # Checks that all errors are small (they definitely should be).
    @test all(vcat(getproperty.(getfield.(wp.wps, :errors), wp.error_estimate)...) .< 10e-9)
    plt = @test_nowarn plot(wp)
    @test length(plt.series_list) == 2

    # Check without appxsol.
    wp = WorkPrecisionSet(static_prob, abstols, reltols, setups;
        names = solnames, numruns = 100, error_estimate = :l2)
    @test all(vcat(getproperty.(getfield.(wp.wps, :errors), wp.error_estimate)...) .< 10e-9)
    plt = @test_nowarn plot(wp)
    @test length(plt.series_list) == 2
end
