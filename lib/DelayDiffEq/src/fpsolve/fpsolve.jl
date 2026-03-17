function OrdinaryDiffEqNonlinearSolve.initial_η(fpsolver::FPSolver, integrator::DDEIntegrator)
    return fpsolver.ηold
end

OrdinaryDiffEqCore.apply_step!(fpsolver::FPSolver, integrator::DDEIntegrator) = nothing

function SciMLBase.postamble!(fpsolver::FPSolver, integrator::DDEIntegrator)
    integrator.stats.nfpiter += fpsolver.iter

    if OrdinaryDiffEqNonlinearSolve.nlsolvefail(fpsolver)
        integrator.stats.nfpconvfail += 1
        @SciMLMessage(
            lazy"Fixed-point iteration failed to converge at t = $(integrator.t) after $(fpsolver.iter) iterations",
            integrator.opts.verbose, :residual_control
        )
    else
        @SciMLMessage(
            lazy"Fixed-point iteration converged at t = $(integrator.t) in $(fpsolver.iter) iterations",
            integrator.opts.verbose, :residual_control
        )
    end
    integrator.force_stepfail = OrdinaryDiffEqNonlinearSolve.nlsolvefail(fpsolver) ||
        integrator.force_stepfail

    return nothing
end
