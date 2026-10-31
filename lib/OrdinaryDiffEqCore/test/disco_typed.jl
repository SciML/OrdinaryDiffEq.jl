using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, SciMLBase, Test
using OrdinaryDiffEqCore: DEVerbosity, SciMLLogging   # SciMLLogging is not a test dependency of its own

# `find_discontinuity` reads the step state through the zero function of a discontinuity problem,
# which `__init` keeps behind a `FunctionWrapper`. For a tuple of callbacks that function's type is
# known from the integrator, so the reads are static and the bracket is a `Tuple`; the root is the
# same as through the wrapper. A `FullSpecialize` problem keeps its callbacks in a tuple; the other
# specializations store them type-erased on Julia 1.12 and later.
f!(du, u, p, t) = (du[1] = 1.0; nothing)
cond(u, t, integrator) = t - 1.0
vcond(out, u, t, integrator) = (out[1] = t - 1.0; out[2] = u[1] - 10.0; nothing)
noop!(integrator) = nothing
vnoop!(integrator, idx) = nothing

# the state `find_discontinuity` sees after a rejected step of size `dt` from `integ.t` (u' = 1)
function straddle!(integ, dt)
    integ.dt = dt
    integ.u .= integ.uprev .+ dt
    return integ
end

const DYNAMIC = (
    :getproperty, :setproperty!, :getindex, :setindex!, :_load_disco_state!,
    :_ode_addsteps!, :solve, :_typevar,
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
            # the stored bracket is left alone; a tuple bracket is built for each solve
            @test dp.tspan == [0.0, 1.0]
            @test OrdinaryDiffEqCore._disco_bracket(dp, 0.25).tspan === (0.0, 0.25)
            # the zero function comes back as the stored object, with its concrete type
            z = OrdinaryDiffEqCore._disco_zero(dp, integ, only(integ.opts.callback.continuous_callbacks))
            @test z === dp.f.f.obj.x
            @test isconcretetype(typeof(z))
            @test only(
                Base.return_types(
                    OrdinaryDiffEqCore._disco_zero,
                    (typeof(dp), typeof(integ), typeof(only(integ.opts.callback.continuous_callbacks)))
                )
            ) === typeof(z)
            # and no dynamic call is left in find_discontinuity for this integrator type
            ci = only(code_typed(OrdinaryDiffEqCore.find_discontinuity, (typeof(integ),); optimize = true))
            dynamic = filter(ci.first.code) do st
                st isa Expr && st.head === :call && st.args[1] isa GlobalRef && st.args[1].name in DYNAMIC
            end
            @test isempty(dynamic)
        end
    end
    # callbacks held in a vector keep the untyped read and find the same root
    vinteg = init(
        prob, Tsit5(); callback = CallbackSet(Any[ContinuousCallback(cond, noop!)], Any[]), dt = 0.5,
        adaptive = false, verbose = DEVerbosity(SciMLLogging.None())
    )
    straddle!(vinteg, 1.5)
    @test OrdinaryDiffEqCore.find_discontinuity(vinteg) ≈ 2 / 3 atol = 1.0e-8
    # end to end: OrdinaryDiffEq's kinked test problem with the event on the kink. With the detection on,
    # a rejected step straddling the event runs `find_discontinuity`, which evaluates the condition
    # through the zero function (its third argument); without it the condition never sees one.
    kink!(du, u, p, t) = (du[1] = u[1] < 1 ? 2u[1] : 5 - 3u[1]; nothing)
    kprob = ODEProblem{true, SciMLBase.FullSpecialize}(kink!, [0.1], (0.0, 1.5))
    zero_func_calls = Ref(0)
    function kcond(u, t, integ)
        integ isa OrdinaryDiffEqCore.zero_func_struct && (zero_func_calls[] += 1)
        return u[1] - 1
    end
    kcb = ContinuousCallback(kcond, noop!)
    # a first step of 1.4 spans the kink (crossed at t ≈ 1.15), so it straddles the event and is rejected
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
    # the kinked problem's solution: u = 0.1exp(2t) up to the kink at τ = log(10)/2 where u = 1,
    # then u' = 5 - 3u from u(τ) = 1
    τ = log(10) / 2
    uref(t) = t <= τ ? 0.1exp(2t) : 5 / 3 - (2 / 3) * exp(-3(t - τ))
    @test isapprox(disco(1.5)[1], uref(1.5); rtol = 1.0e-5, atol = 1.0e-8)
    @test isapprox(disco.u[end][1], uref(1.5); rtol = 1.0e-5, atol = 1.0e-8)
end

@testset "several callbacks, repeated detection" begin
    # time conditions (`straddle!` sets `u` but not the interpolant); t - 0.8 and t - 0.3 both cross on [0, 1]
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
        # a later, shorter step: the root at t = 0.3 is θ = 0.6 of a step of 0.5, whatever the last bracket was
        straddle!(integ, 0.5)
        @test OrdinaryDiffEqCore.find_discontinuity(integ) ≈ 0.6 atol = 1.0e-8
        @test all(dp -> dp.tspan == [0.0, 1.0], dps)
    end
end
