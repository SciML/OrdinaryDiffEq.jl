module OrdinaryDiffEqRosenbrockReverseDiffExt

using OrdinaryDiffEqRosenbrock, ReverseDiff
using ADTypes: AutoFiniteDiff

function OrdinaryDiffEqRosenbrock.check_reversediff_compatibility(
        u::AbstractArray{<:ReverseDiff.TrackedReal}, alg
    )
    if alg.autodiff isa AutoFiniteDiff && alg.linsolve === nothing
        throw(ArgumentError("ReverseDiff through an in-place Rosenbrock solve with AutoFiniteDiff and the default linear solver is unsupported. Use ForwardDiff or a SciMLSensitivity sensitivity algorithm instead of SensitivityADPassThrough()."))
    end
    return nothing
end

end
