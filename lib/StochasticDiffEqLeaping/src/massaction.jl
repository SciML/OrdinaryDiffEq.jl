struct MassActionRates{M}
    jump::M
end

(rate::MassActionRates)(out, u, p, t) =
    JumpProcesses.massaction_rates!(out, rate.jump, u)
function (rate::MassActionRates)(u, p, t)
    out = similar(rate.jump.scaled_rates, promote_type(eltype(u), eltype(rate.jump.scaled_rates)))
    return JumpProcesses.massaction_rates!(out, rate.jump, u)
end

function StochasticDiffEqCore.jump_noise_data(
        alg::Union{TauLeaping, CaoTauLeaping, ImplicitTauLeaping, ThetaTrapezoidalTauLeaping},
        prob, u, p, t
    )
    if !(prob isa JumpProcesses.JumpProblem) ||
            !(prob.aggregator isa JumpProcesses.PureLeaping) ||
            JumpProcesses.get_num_majumps(prob.massaction_jump) == 0
        return invoke(
            StochasticDiffEqCore.jump_noise_data,
            Tuple{Any, Any, Any, Any, Any}, alg, prob, u, p, t
        )
    end
    prob.prob isa SciMLBase.DiscreteProblem ||
        throw(ArgumentError("Mass-action leaping requires a DiscreteProblem with PureLeaping()"))
    prob.regular_jump === nothing || throw(
        ArgumentError(
            "Combining MassActionJump and RegularJump is not supported by this solver"
        )
    )
    isempty(prob.constant_jumps) && isempty(prob.variable_jumps) || throw(
        ArgumentError(
            "Mass-action leaping requires a pure mass-action problem"
        )
    )
    jump = prob.massaction_jump
    rate_type = promote_type(float(eltype(u)), float(eltype(jump.scaled_rates)))
    return (;
        jump_prototype = zeros(rate_type, JumpProcesses.get_num_majumps(jump)), c = jump,
        rate_constants = copy(jump.scaled_rates), rate = MassActionRates(jump),
        iip = SciMLBase.isinplace(prob),
    )
end

leaping_change(c, args...) = c(args...)
function leaping_change(jump::JumpProcesses.MassActionJump, du, u, p, t, counts, mark)
    return JumpProcesses.massaction_stoichiometry_mul!(du, jump, counts)
end
function leaping_change(jump::JumpProcesses.MassActionJump, u, p, t, counts, mark)
    du = similar(u)
    JumpProcesses.massaction_stoichiometry_mul!(du, jump, counts)
    return convert(typeof(u), du)
end

struct MassActionDrift{M, IIP}
    jump::M
end
(drift::MassActionDrift{M, true})(du, u, p, t) where {M} =
    JumpProcesses.massaction_drift!(du, drift.jump, u)
function (drift::MassActionDrift{M, false})(u, p, t) where {M}
    du = similar(u)
    JumpProcesses.massaction_drift!(du, drift.jump, u)
    return convert(typeof(u), du)
end
