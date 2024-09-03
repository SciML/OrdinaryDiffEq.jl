@doc generic_solver_docstring(
    "Feagin's 10th-order method.", "Feagin10", "Explicit Runge-Kutta Method. ",
    """@article{feagin2012high,
    title={High-order explicit Runge-Kutta methods using m-symmetry},
    author={Feagin, Terry},
    year={2012},
    publisher={Neural, Parallel \\& Scientific Computations}
    }""",
    """
     - `stage_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    """
    step_limiter! = OrdinaryDiffEq.trivial_limiter!,
    """
)
Base.@kwdef struct Feagin10{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
    step_limiter!::StepLimiter = trivial_limiter!
end

@doc generic_solver_docstring(
    "Feagin's 12th-order method.", "Feagin12", "Explicit Runge-Kutta Method. ",
    """@article{feagin2012high,
    title={High-order explicit Runge-Kutta methods using m-symmetry},
    author={Feagin, Terry},
    year={2012},
    publisher={Neural, Parallel \\& Scientific Computations}
    }""",
    """
     - `stage_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    """
    step_limiter! = OrdinaryDiffEq.trivial_limiter!,
    """
)
Base.@kwdef struct Feagin12{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
    step_limiter!::StepLimiter = trivial_limiter!
end

@doc generic_solver_docstring(
    "Feagin's 14th-order method.", "Feagin14", "Explicit Runge-Kutta Method. ",
    """@article{feagin2009explicit,
    title={An Explicit Runge-Kutta Method of Order Fourteen},
    author={Feagin, Terry},
    year={2009},
    publisher={Numerical Algorithms}
    }""",
    """
     - `stage_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    """
    step_limiter! = OrdinaryDiffEq.trivial_limiter!,
    """
)
Base.@kwdef struct Feagin14{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
    step_limiter!::StepLimiter = trivial_limiter!
end
