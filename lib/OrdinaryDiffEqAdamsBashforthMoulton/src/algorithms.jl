# Adams Bashforth and Adams moulton methods

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

AB3: Adams-Bashforth Explicit Method
The 3-step third order multistep method. Ralston's Second Order Method is used to calculate starting values.
"""
struct AB3 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

AB4: Adams-Bashforth Explicit Method
The 4-step fourth order multistep method. Runge-Kutta method of order 4 is used to calculate starting values.
"""
struct AB4 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

AB5: Adams-Bashforth Explicit Method
The 3-step third order multistep method. Ralston's Second Order Method is used to calculate starting values.
"""
struct AB5 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

ABM32: Adams-Bashforth Explicit Method
It is third order method. In ABM32, AB3 works as predictor and Adams Moulton 2-steps method works as Corrector.
Ralston's Second Order Method is used to calculate starting values.
"""
struct ABM32 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

ABM43: Adams-Bashforth Explicit Method
It is fourth order method. In ABM43, AB4 works as predictor and Adams Moulton 3-steps method works as Corrector.
Runge-Kutta method of order 4 is used to calculate starting values.
"""
struct ABM43 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

ABM54: Adams-Bashforth Explicit Method
It is fifth order method. In ABM54, AB5 works as predictor and Adams Moulton 4-steps method works as Corrector.
Runge-Kutta method of order 4 is used to calculate starting values.
"""
struct ABM54 <: OrdinaryDiffEqAlgorithm end

# Variable Step Size Adams methods

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCAB3: Adaptive step size Adams explicit Method
The 3rd order Adams method. Bogacki-Shampine 3/2 method is used to calculate starting values.
"""
struct VCAB3 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCAB4: Adaptive step size Adams explicit Method
The 4th order Adams method. Runge-Kutta 4 is used to calculate starting values.
"""
struct VCAB4 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCAB5: Adaptive step size Adams explicit Method
The 5th order Adams method. Runge-Kutta 4 is used to calculate starting values.
"""
struct VCAB5 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCABM3: Adaptive step size Adams explicit Method
The 3rd order Adams-Moulton method. Bogacki-Shampine 3/2 method is used to calculate starting values.
"""
struct VCABM3 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCABM4: Adaptive step size Adams explicit Method
The 4th order Adams-Moulton method. Runge-Kutta 4 is used to calculate starting values.
"""
struct VCABM4 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCABM5: Adaptive step size Adams explicit Method
The 5th order Adams-Moulton method. Runge-Kutta 4 is used to calculate starting values.
"""
struct VCABM5 <: OrdinaryDiffEqAdaptiveAlgorithm end

# Variable Order and Variable Step Size Adams methods

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCABM: Adaptive step size Adams explicit Method
An adaptive order adaptive time Adams Moulton method.
It uses an order adaptivity algorithm is derived from Shampine's DDEABM.
"""
struct VCABM <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm end

const MultistepAlgorithms = Union{
    AB3, AB4, AB5, ABM32, ABM43, ABM54}