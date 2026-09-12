using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, Test
using OrdinaryDiffEqCore: fix_dt_at_bounds!
using OrdinaryDiffEqCore.DiffEqBase: timedepentdtmin

# `fix_dt_at_bounds!` clamps |dt| into [dtmin, dtmax]. `timedepentdtmin` returns a
# positive magnitude, so the lower clamp for a backward solve has to compare
# against `-dtmin`; `min(integrator.dt, dtmin)` is a no-op whenever dt < 0 < dtmin,
# which silently removed the floor for every tdir < 0 solve.

@testset "dtmin direction handling" begin
@testset "dtmin floor applies in both directions" begin
    for (t0, t1) in ((0.0, 1.0), (1.0, 0.0))
        prob = ODEProblem((u, p, t) -> -u, 1.0, (t0, t1))
        integ = init(prob, Tsit5(); dtmin = 1.0e-3, force_dtmin = true)
        dtmin = timedepentdtmin(integ)
        integ.dt = integ.tdir * dtmin / 100 # far below the floor
        fix_dt_at_bounds!(integ)
        @test abs(integ.dt) >= dtmin
        @test signbit(integ.dt) == signbit(integ.tdir)
    end
end

# End-to-end consequence: on a tspan symmetric about zero the mirrored problem
# `g(u, p, t) = -f(u, p, -t)` traverses the same trajectory backwards, and because
# IEEE negation and round-to-nearest are sign-symmetric the two solves agree
# bitwise. A floor applied in only one direction breaks that.
@testset "time-reversal symmetry at the dtmin floor" begin
    f = (u, p, t) -> -50.0 * u
    g = (u, p, t) -> 50.0 * u # -f(u, p, -t) for this autonomous f
    kw = (
        abstol = 1.0e-12, reltol = 1.0e-12, dtmin = 1.0e-2,
        force_dtmin = true, save_everystep = false,
    )
    dts(prob) = (i = init(prob, Tsit5(); kw...); [abs(Float64(i.dt)) for _ in i])
    fwd = dts(ODEProblem(f, 1.0, (-1.0, 1.0)))
    bwd = dts(ODEProblem(g, 1.0, (1.0, -1.0)))
    @test !isempty(fwd)
    @test fwd == bwd
end
end
