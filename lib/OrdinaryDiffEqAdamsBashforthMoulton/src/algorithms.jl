# Adams Bashforth and Adams moulton methods
reference = """E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
            Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
            https://doi.org/10.1007/978-3-540-78862-1"""

@doc generic_solver_docstring("The 3-step third order multistep method.
        Ralston's Second Order Method is used to calculate starting values.",
    "AB3",
    "Adams-Bashforth Explicit Method",
    reference,
    "",
    "")
struct AB3{Thread} <: OrdinaryDiffEqAlgorithm
    thread::Thread
end

Base.@kwdef struct AB3{Thread} <: OrdinaryDiffEqAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("The 4-step fourth order multistep method.
    Runge-Kutta method of order 4 is used to calculate starting values.",
    "AB4",
    "Adams-Bashforth Explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct AB4{Thread} <: OrdinaryDiffEqAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("The 5-step fifth order multistep method.
    Ralston's 3rd order Runge-Kutta method is used to calculate starting values.",
    "AB5",
    "Adams-Bashforth Explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct AB5{Thread} <: OrdinaryDiffEqAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("It is third order method.
    In ABM32, AB3 works as predictor and Adams Moulton 2-steps method works as Corrector.
    Ralston's Second Order Method is used to calculate starting values.",
    "ABM32",
    "Adams-Bashforth Explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct ABM32{Thread} <: OrdinaryDiffEqAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("It is fourth order method.
    In ABM43, AB4 works as predictor and Adams Moulton 3-steps method works as Corrector.
    Runge-Kutta method of order 4 is used to calculate starting values.",
    "ABM43",
    "Adams-Bashforth Explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct ABM43{Thread} <: OrdinaryDiffEqAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("It is fifth order method.
    In ABM54, AB5 works as predictor and Adams Moulton 4-steps method works as Corrector.
    Runge-Kutta method of order 4 is used to calculate starting values.",
    "ABM54",
    "Adams-Bashforth Explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct ABM54{Thread} <: OrdinaryDiffEqAlgorithm
    thread::Thread = False()
end

# Variable Step Size Adams methods

@doc generic_solver_docstring("The 3rd order Adams method.
    Bogacki-Shampine 3/2 method is used to calculate starting values.",
    "VCAB3",
    "Adams explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct VCAB3{Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("The 4th order Adams method.
    Runge-Kutta 4 is used to calculate starting values.",
    "VCAB4",
    "Adams explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct VCAB4{Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("The 5th order Adams method.
    Runge-Kutta 4 is used to calculate starting values.",
    "VCAB5",
    "Adams explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct VCAB5{Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("The 3rd order Adams-Moulton method.
    Bogacki-Shampine 3/2 method is used to calculate starting values.",
    "VCABM3",
    "Adams explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct VCABM3{Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("The 4th order Adams-Moulton method.
    Runge-Kutta 4 is used to calculate starting values.",
    "VCABM4",
    "Adams explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct VCABM4{Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    thread::Thread = False()
end

@doc generic_solver_docstring("The 5th order Adams-Moulton method.
    Runge-Kutta 4 is used to calculate starting values.",
    "VCABM5",
    "Adams explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct VCABM5{Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    thread::Thread = False()
end

# Variable Order and Variable Step Size Adams methods

@doc generic_solver_docstring("An adaptive order adaptive time Adams Moulton method.
    It uses an order adaptivity algorithm is derived from Shampine's DDEABM.",
    "VCABM",
    "adaptive order Adams explicit Method",
    reference,
    "",
    "")
Base.@kwdef struct VCABM{Thread} <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm
    thread::Thread = False()
end
