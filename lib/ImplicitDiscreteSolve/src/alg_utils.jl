function SciMLBase.isautodifferentiable(alg::IDSolve)
    true
end
function SciMLBase.allows_arbitrary_number_types(alg::IDSolve)
    true
end
function SciMLBase.allowscomplex(alg::IDSolve)
    true
end

SciMLBase.isdiscrete(alg::IDSolve) = true

isfsal(alg::IDSolve) = false
alg_order(alg::IDSolve) = 1
beta2_default(alg::IDSolve) = 0
beta1_default(alg::IDSolve, beta2) = 0

dt_required(alg::IDSolve) = false
isdiscretealg(alg::IDSolve) = true

isadaptive(alg::IDSolve) = true

# @concrete struct ConvergenceRateTracing
#     inner_tracing
# end

# @concrete struct ConvergenceRateTraceTrick
#     incrementL2norms
#     residualL2norms
#     trace_wrapper
# end

# function NonlinearSolveBase.init_nonlinearsolve_trace(
#         prob, alg::IDSolve, u, fu, J, δu;
#         trace_level::ConvergenceRateTracing, kwargs... # This kind of dispatch does not work. Need to figure out a different way.
# )
#     inner_trace = NonlinearSolveBase.init_nonlinearsolve_trace(
#         prob, alg, u, fu, J, δu;
#         trace_level.inner_tracing, kwargs...
#     )

#     return ConvergenceRateTraceTrick(eltype(δu)[], eltype(fu)[], inner_trace)
# end


