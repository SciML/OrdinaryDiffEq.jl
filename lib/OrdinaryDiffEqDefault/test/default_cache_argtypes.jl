using OrdinaryDiffEqDefault, OrdinaryDiffEqCore, ADTypes, Test

f(du, u, p, t) = (du .= -u)
prob = ODEProblem(f, [1.0], (0.0, 1.0))
integ = init(prob, DefaultODEAlgorithm(autodiff = AutoForwardDiff(chunksize = 1)))

@test isconcretetype(fieldtype(typeof(integ.cache), :args))

ci, _ = only(
    code_typed(
        OrdinaryDiffEqCore.init_ith_default_cache,
        (typeof(integ.cache), typeof(integ.alg.algs), Int); optimize = true
    )
)
dynamic_alg_cache = filter(ci.code) do st
    st isa Expr && st.head === :call && any(a -> a isa GlobalRef && a.name === :alg_cache, st.args)
end
@test isempty(dynamic_alg_cache)

sol = solve!(integ)
@test sol.retcode == ReturnCode.Success
@test sol.u[end][1] ≈ exp(-1) rtol = 1.0e-3
