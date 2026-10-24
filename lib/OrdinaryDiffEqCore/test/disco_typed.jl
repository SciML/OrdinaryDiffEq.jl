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
    # end to end with the detection on: the time root is hit and the solve succeeds
    sols = map((SciMLLogging.Standard(), SciMLLogging.None())) do v
        integ = init(
            prob, Tsit5(); callback = ContinuousCallback(cond, noop!), verbose = DEVerbosity(v),
            controller = OrdinaryDiffEqCore.PIController(Float64, Tsit5(); discontinuity_detection = true)
        )
        solve!(integ)
    end
    for sol in sols
        @test sol.retcode == ReturnCode.Success
        @test any(t -> abs(t - 1.0) < 1.0e-6, sol.t)
    end
    @test sols[1].t == sols[2].t && sols[1].u == sols[2].u
end
