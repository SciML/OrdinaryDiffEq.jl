using OrdinaryDiffEqDefault, OrdinaryDiffEqCore, ADTypes, ForwardDiff, StaticArrays, Test

f(du, u, p, t) = (du .= -u)
prob = ODEProblem(f, [1.0], (0.0, 1.0))
integ = init(prob, DefaultODEAlgorithm(autodiff = AutoForwardDiff(chunksize = 1)))

@test isconcretetype(fieldtype(typeof(integ.cache), :args))
# the element-type arguments are stored as singletons, so their types alone determine `alg_cache`'s method
@test all(Base.issingletontype, fieldtype(typeof(integ.cache), :args).parameters[3:5])

function dynamic_alg_cache_calls(cache, algs)
    ci, _ = only(
        code_typed(
            OrdinaryDiffEqCore.init_ith_default_cache,
            (typeof(cache), typeof(algs), Int); optimize = true
        )
    )
    return filter(ci.code) do st
        st isa Expr && st.head === :call && any(a -> a isa GlobalRef && a.name === :alg_cache, st.args)
    end
end
@test isempty(dynamic_alg_cache_calls(integ.cache, integ.alg.algs))

sol = solve!(integ)
@test sol.retcode == ReturnCode.Success
@test sol.u[end][1] ≈ exp(-1) rtol = 1.0e-3

# Every member of the cache, built lazily or eagerly, for in-place, out-of-place, static and dual states.
# The element types handed to `alg_cache` are those of the state and time, and building a member a second
# time keeps the first one.
f_oop(u, p, t) = -p * u
decay_problems = (
    "in-place" => ODEProblem(f, [1.0], (0.0, 1.0)),
    "out-of-place" => ODEProblem(f_oop, [1.0], (0.0, 1.0), 1.0),
    "static" => ODEProblem(f_oop, SA[1.0], (0.0, 1.0), 1.0),
    "dual" => ODEProblem(f, [ForwardDiff.Dual{typeof(ForwardDiff.Tag(f, Float64))}(1.0, 1.0)], (0.0, 1.0)),
)
@testset "DefaultCache members: $name" for (name, prob) in decay_problems
    alg = DefaultODEAlgorithm(autodiff = AutoForwardDiff(chunksize = 1))
    integ = init(prob, alg; abstol = 1.0e-10, reltol = 1.0e-10)
    cache = integ.cache
    @test isempty(dynamic_alg_cache_calls(cache, integ.alg.algs))
    args = OrdinaryDiffEqCore.default_cache_args(cache.args...)
    @test args[3] === eltype(prob.u0)
    @test args[4] === eltype(prob.u0)
    @test args[5] === Float64
    for k in 1:6
        member = Symbol(:cache, k)
        OrdinaryDiffEqCore.init_ith_default_cache(cache, integ.alg.algs, k)
        @test isdefined(cache, member)
        built = getfield(cache, member)
        @test built isa fieldtype(typeof(cache), member)
        OrdinaryDiffEqCore.init_ith_default_cache(cache, integ.alg.algs, k)
        @test getfield(cache, member) === built
    end
    sol = solve!(integ)
    @test sol.retcode == ReturnCode.Success
    @test ForwardDiff.value(sol.u[end][1]) ≈ exp(-1) rtol = 1.0e-6
    if name == "dual"
        @test ForwardDiff.partials(sol.u[end][1])[1] ≈ exp(-1) rtol = 1.0e-6
    end
end

# d/dp of exp(-p) at p = 1, differentiated through the default algorithm's cache construction
@test ForwardDiff.derivative(
    p -> solve(
        ODEProblem(f_oop, 1.0, (0.0, 1.0), p), DefaultODEAlgorithm();
        abstol = 1.0e-10, reltol = 1.0e-10
    ).u[end],
    1.0
) ≈ -exp(-1) rtol = 1.0e-6
