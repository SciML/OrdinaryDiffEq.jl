"""
    AMF(AlgType; kwargs...)

Approximate Matrix Factorization wrapper for Rosenbrock-W methods.

Wraps any Rosenbrock-W algorithm type (e.g. `ROS34PW1a`, `Rosenbrock23`) and
automatically uses `SciMLOpFactorization()` as the linear solver, so that the
`W_prototype` SciMLOperator product is factorized per-factor instead of being
concretized into a dense matrix.

The ODEFunction must be built with `build_amf_function` (or manually with
`jac_prototype` and `W_prototype` SciMLOperators).

# Usage

```julia
func = build_amf_function(f!; n = 20, jac_upper = fjac_upper, jac_lower = fjac_lower)
prob = ODEProblem(func, u0, tspan)
sol = solve(prob, AMF(ROS34PW1a))
```

Additional keyword arguments are forwarded to the inner algorithm constructor:

```julia
sol = solve(prob, AMF(ROS34PW1a; autodiff = AutoFiniteDiff()))
```
"""
struct AMF{AlgType, K} <: OrdinaryDiffEqAlgorithm
    algtype::Type{AlgType}
    kwargs::K
end

AMF(AlgType::Type; kwargs...) = AMF(AlgType, kwargs)

function _build_inner_alg(amf::AMF)
    return amf.algtype(; linsolve = SciMLOpFactorization(), amf.kwargs...)
end

function SciMLBase.__solve(prob::SciMLBase.AbstractODEProblem, amf::AMF, args...; kwargs...)
    return SciMLBase.__solve(prob, _build_inner_alg(amf), args...; kwargs...)
end

function SciMLBase.__init(prob::SciMLBase.AbstractODEProblem, amf::AMF, args...; kwargs...)
    return SciMLBase.__init(prob, _build_inner_alg(amf), args...; kwargs...)
end
