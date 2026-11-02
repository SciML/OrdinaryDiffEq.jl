using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, SciMLBase, StaticArrays, Test
using OrdinaryDiffEqCore: DEVerbosity, SciMLLogging

# FullSpecialize keeps callbacks in a tuple (typed zero-function reads); a vector of callbacks is read untyped.
f!(du, u, p, t) = (du[1] = 1.0; nothing)
cond(u, t, integrator) = t - 1.0
vcond(out, u, t, integrator) = (out[1] = t - 1.0; out[2] = u[1] - 10.0; nothing)
noop!(integrator) = nothing
vnoop!(integrator, idx) = nothing

# the state after a rejected step of size `dt` from `integ.t` (u' = 1)
function straddle!(integ, dt)
    integ.dt = dt
    integ.u .= integ.uprev .+ dt
    return integ
end

const DYNAMIC = (
    :getproperty, :setproperty!, :getindex, :setindex!, :_load_disco_state!,
    :_set_disco_outs!, :_ode_addsteps!, :solve, :_typevar,
)

@testset "typed discontinuity detection" begin
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(f!, [0.0], (0.0, 2.0))
    for (cb, name) in (
            (ContinuousCallback(cond, noop!), "scalar"),
            (VectorContinuousCallback(vcond, vnoop!, 2), "vector"),
        )
        @testset "$name" begin
            integ = init(
                prob, Tsit5(); callback = cb, dt = 0.5, adaptive = false,
                verbose = DEVerbosity(SciMLLogging.None())
            )
            straddle!(integ, 1.5)
            θ = OrdinaryDiffEqCore.find_discontinuity(integ)
            @test θ ≈ 2 / 3 atol = 1.0e-8              # the root of t - 1 on the step [0, 1.5]
            @test integ.is_disco_step
            dp = only(OrdinaryDiffEqCore.get_disco_probs(integ.controller_cache))
            @test dp.tspan == [0.0, 1.0]
            @test OrdinaryDiffEqCore._disco_bracket(dp, 0.25).tspan === (0.0, 0.25)
            z = OrdinaryDiffEqCore._disco_zero(dp, integ, only(integ.opts.callback.continuous_callbacks))
            @test z === dp.f.f.obj.x
            @test isconcretetype(typeof(z))
            @test only(
                Base.return_types(
                    OrdinaryDiffEqCore._disco_zero,
                    (typeof(dp), typeof(integ), typeof(only(integ.opts.callback.continuous_callbacks)))
                )
            ) === typeof(z)
            ci = only(code_typed(OrdinaryDiffEqCore.find_discontinuity, (typeof(integ),); optimize = true))
            dynamic = filter(ci.first.code) do st
                st isa Expr && st.head === :call && st.args[1] isa GlobalRef && st.args[1].name in DYNAMIC
            end
            @test isempty(dynamic)
        end
    end
    vinteg = init(
        prob, Tsit5(); callback = CallbackSet(Any[ContinuousCallback(cond, noop!)], Any[]), dt = 0.5,
        adaptive = false, verbose = DEVerbosity(SciMLLogging.None())
    )
    straddle!(vinteg, 1.5)
    @test OrdinaryDiffEqCore.find_discontinuity(vinteg) ≈ 2 / 3 atol = 1.0e-8
    kink!(du, u, p, t) = (du[1] = u[1] < 1 ? 2u[1] : 5 - 3u[1]; nothing)
    kprob = ODEProblem{true, SciMLBase.FullSpecialize}(kink!, [0.1], (0.0, 1.5))
    zero_func_calls = Ref(0)
    function kcond(u, t, integ)
        integ isa OrdinaryDiffEqCore.zero_func_struct && (zero_func_calls[] += 1)
        return u[1] - 1
    end
    kcb = ContinuousCallback(kcond, noop!)
    # a first step of 1.4 straddles the kink at t ≈ 1.15 and is rejected
    disco = solve(
        kprob, Tsit5(); callback = kcb, reltol = 1.0e-6, dt = 1.4,
        controller = OrdinaryDiffEqCore.PIController(Tsit5(); discontinuity_detection = true)
    )
    @test disco.retcode == ReturnCode.Success
    @test zero_func_calls[] > 0
    zero_func_calls[] = 0
    plain = solve(kprob, Tsit5(); callback = kcb, reltol = 1.0e-6, dt = 1.4)
    @test zero_func_calls[] == 0
    @test disco.stats.nreject <= plain.stats.nreject
    # analytic: u = 0.1exp(2t) up to τ = log(10)/2 where u = 1, then u' = 5 - 3u
    τ = log(10) / 2
    uref(t) = t <= τ ? 0.1exp(2t) : 5 / 3 - (2 / 3) * exp(-3(t - τ))
    @test isapprox(disco(1.5)[1], uref(1.5); rtol = 1.0e-5, atol = 1.0e-8)
    @test isapprox(disco.u[end][1], uref(1.5); rtol = 1.0e-5, atol = 1.0e-8)
end

@testset "several callbacks, repeated detection" begin
    # time conditions: `straddle!` sets `u` but not the interpolant
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(f!, [0.0], (0.0, 2.0))
    c08(u, t, integrator) = t - 0.8
    c03(u, t, integrator) = t - 0.3
    for cbs in ((c08, c03),)
        integ = init(
            prob, Tsit5(); callback = CallbackSet(map(c -> ContinuousCallback(c, noop!), cbs)...),
            dt = 0.5, adaptive = false, verbose = DEVerbosity(SciMLLogging.None())
        )
        @test integ.opts.callback.continuous_callbacks isa Tuple
        straddle!(integ, 1.0)
        @test OrdinaryDiffEqCore.find_discontinuity(integ) ≈ 0.3 atol = 1.0e-8
        dps = OrdinaryDiffEqCore.get_disco_probs(integ.controller_cache)
        @test length(dps) == 2
        @test all(dp -> dp.tspan == [0.0, 1.0], dps)
        straddle!(integ, 0.5)
        @test OrdinaryDiffEqCore.find_discontinuity(integ) ≈ 0.6 atol = 1.0e-8
        @test all(dp -> dp.tspan == [0.0, 1.0], dps)
    end
end

warm_bytes(integ) = (OrdinaryDiffEqCore.find_discontinuity(integ); @allocated OrdinaryDiffEqCore.find_discontinuity(integ))

@testset "scalar callback end values" begin
    cb = ContinuousCallback(cond, noop!)
    @test only(
        Base.return_types(OrdinaryDiffEqCore._disco_outs, (Vector{Float64}, typeof(cb)))
    ) === Tuple{Nothing, Float64, Float64}
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(f!, [0.0], (0.0, 2.0))
    integ = init(
        prob, Tsit5(); callback = cb, dt = 0.5, adaptive = false,
        verbose = DEVerbosity(SciMLLogging.None())
    )
    straddle!(integ, 1.5)
    @test OrdinaryDiffEqCore.find_discontinuity(integ) ≈ 2 / 3 atol = 1.0e-8
    z = only(OrdinaryDiffEqCore.get_disco_probs(integ.controller_cache)).f.f.obj.x
    @test fieldtype(typeof(z), :out_low) === Float64
    @test fieldtype(typeof(z), :out_high) === Float64
    @test z(0.0, SciMLBase.NullParameters()) === cond(integ.uprev, 0.0, integ) === -1.0
    @test z(1.0, SciMLBase.NullParameters()) === cond(integ.u, 1.5, integ) === 0.5
    @test warm_bytes(integ) == 0
end

# u' = (1, 2) from the origin: the first component crosses 1 at θ = 2/3 of a step of 1.5.
# The conditions take the state type of the problem; the step is a real one, so the
# zero function interpolates through its stages.
static_f(u, p, t) = SVector(1.0, 2.0)
static_cond(u::SVector, t, integrator) = u[1] - 1.0
function static_vcond(out, u::SVector, t, integrator)
    out[1] = u[1] - 1.0
    out[2] = u[2] - 10.0
    return nothing
end

@testset "static array state, $name, $spec" for spec in (SciMLBase.FullSpecialize, SciMLBase.AutoSpecialize),
        (name, cb) in (
            ("scalar", ContinuousCallback(static_cond, noop!)),
            ("vector", VectorContinuousCallback(static_vcond, vnoop!, 2)),
        )

    prob = ODEProblem{false, spec}(static_f, SVector(0.0, 0.0), (0.0, 2.0))
    integ = init(
        prob, Tsit5(); callback = cb, dt = 1.5, adaptive = false,
        verbose = DEVerbosity(SciMLLogging.None())
    )
    OrdinaryDiffEqCore.perform_step!(integ, integ.cache)
    @test OrdinaryDiffEqCore.find_discontinuity(integ) ≈ 2 / 3 atol = 1.0e-8
    @test integ.is_disco_step
    dp = only(OrdinaryDiffEqCore.get_disco_probs(integ.controller_cache))
    @test dp.f.f.obj.x.u₁ isa SVector{2, Float64}
    if spec === SciMLBase.FullSpecialize
        icb = only(integ.opts.callback.continuous_callbacks)
        @test only(
            Base.return_types(
                OrdinaryDiffEqCore._disco_zero, (typeof(dp), typeof(integ), typeof(icb))
            )
        ) === typeof(dp.f.f.obj.x)
        @test warm_bytes(integ) == 0
    end

    kink(u, p, t) = SVector(u[1] < 1 ? 2u[1] : 5 - 3u[1], 0.0)
    τ = log(10) / 2
    kprob = ODEProblem{false, spec}(kink, SVector(0.1, 0.0), (0.0, 1.5))
    disco = solve(
        kprob, Tsit5(); callback = cb, reltol = 1.0e-6, dt = 1.4,
        controller = OrdinaryDiffEqCore.PIController(Tsit5(); discontinuity_detection = true)
    )
    @test disco.retcode == ReturnCode.Success
    @test isapprox(disco.u[end][1], 5 / 3 - (2 / 3) * exp(-3(1.5 - τ)); rtol = 1.0e-5)
end

@testset "mutable static array state keeps the interpolation buffer" begin
    mutable_f!(du, u, p, t) = (du[1] = 1.0; du[2] = 2.0; nothing)
    mutable_cond(u::MVector, t, integrator) = u[1] - 1.0
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(mutable_f!, MVector(0.0, 0.0), (0.0, 2.0))
    integ = init(
        prob, Tsit5(); callback = ContinuousCallback(mutable_cond, noop!), dt = 1.5,
        adaptive = false, verbose = DEVerbosity(SciMLLogging.None())
    )
    OrdinaryDiffEqCore.perform_step!(integ, integ.cache)
    @test OrdinaryDiffEqCore.find_discontinuity(integ) ≈ 2 / 3 atol = 1.0e-8
    z = only(OrdinaryDiffEqCore.get_disco_probs(integ.controller_cache)).f.f.obj.x
    @test z.u₁ isa MVector{2, Float64}
    @test z.u₁ !== integ.u
end
