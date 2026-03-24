@static if isdefined(SciMLBase, :DEStats)
    const Stats = SciMLBase.DEStats
else
    """
    $(TYPEDEF)

    Statistics from the differential equation solver about the solution process.

    ## Fields

    - nf: Number of function evaluations. If the differential equation is a split function,
      such as a `SplitFunction` for implicit-explicit (IMEX) integration, then `nf` is the
      number of function evaluations for the first function (the implicit function)
    - nf2: If the differential equation is a split function, such as a `SplitFunction`
      for implicit-explicit (IMEX) integration, then `nf2` is the number of function
      evaluations for the second function, i.e. the function treated explicitly. Otherwise
      it is zero.
    - nw: The number of W=I-gamma*J (or W=I/gamma-J) matrices constructed during the solving
      process.
    - nsolve: The number of linear solves `W\b` required for the integration.
    - njacs: Number of Jacobians calculated during the integration.
    - nnonliniter: Total number of iterations for the nonlinear solvers.
    - nnonlinconvfail: Number of nonlinear solver convergence failures.
    - ncondition: Number of calls to the condition function for callbacks.
    - naccept: Number of accepted steps.
    - nreject: Number of rejected steps.
    - maxeig: Maximum eigenvalue over the solution. This is only computed if the
      method is an auto-switching algorithm.
    """
    mutable struct Stats
        nf::Int
        nf2::Int
        nw::Int
        nsolve::Int
        njacs::Int
        nnonliniter::Int
        nnonlinconvfail::Int
        ncondition::Int
        naccept::Int
        nreject::Int
        maxeig::Float64
    end

    Base.@deprecate_binding DEStats Stats false

    Stats(x::Int = -1) = Stats(x, x, x, x, x, x, x, x, x, x, 0.0)

    function Base.show(io::IO, s::Stats)
        println(io, summary(s))
        @printf io "%-50s %-d\n" "Number of function 1 evaluations:" s.nf
        @printf io "%-50s %-d\n" "Number of function 2 evaluations:" s.nf2
        @printf io "%-50s %-d\n" "Number of W matrix evaluations:" s.nw
        @printf io "%-50s %-d\n" "Number of linear solves:" s.nsolve
        @printf io "%-50s %-d\n" "Number of Jacobians created:" s.njacs
        @printf io "%-50s %-d\n" "Number of nonlinear solver iterations:" s.nnonliniter
        @printf io "%-50s %-d\n" "Number of nonlinear solver convergence failures:" s.nnonlinconvfail
        @printf io "%-50s %-d\n" "Number of rootfind condition calls:" s.ncondition
        @printf io "%-50s %-d\n" "Number of accepted steps:" s.naccept
        @printf io "%-50s %-d" "Number of rejected steps:" s.nreject
        iszero(s.maxeig) || @printf io "\n%-50s %-d" "Maximum eigenvalue recorded:" s.maxeig
    end
end
