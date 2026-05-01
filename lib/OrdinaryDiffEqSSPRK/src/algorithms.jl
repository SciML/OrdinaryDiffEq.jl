@doc explicit_rk_docstring(
    "A third-order, five-stage explicit strong stability preserving (SSP) low-storage method.
Fixed timestep only.",
    "SSPRK53_2N2",
    references = "Higueras and T. Roldán.
    New third order low-storage SSP explicit Runge–Kutta methods
    arXiv:1809.04807v1."
)
Base.@kwdef struct SSPRK53_2N2{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A second-order, two-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK22",
    references = "Shu, Chi-Wang, and Stanley Osher.
    Efficient implementation of essentially non-oscillatory shock-capturing schemes.
    Journal of Computational Physics 77.2 (1988): 439-471.
    https://doi.org/10.1016/0021-9991(88)90177-5"
)
Base.@kwdef struct SSPRK22{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, five-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK53",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207"
)
Base.@kwdef struct SSPRK53{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, six-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK63",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207"
)
Base.@kwdef struct SSPRK63{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, eight-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK83",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207"
)
Base.@kwdef struct SSPRK83{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, four-stage explicit strong stability preserving (SSP) method.",
    "SSPRK43",
    references = """Optimal third-order explicit SSP method with four stages discovered by

      - J. F. B. M. Kraaijevanger.
        "Contractivity of Runge-Kutta methods."
        In: BIT Numerical Mathematics 31.3 (1991), pp. 482–528.
        [DOI: 10.1007/BF01933264](https://doi.org/10.1007/BF01933264).

    Embedded method constructed by

      - Sidafa Conde, Imre Fekete, John N. Shadid.
        "Embedded error estimation and adaptive step-size control for
        optimal explicit strong stability preserving Runge–Kutta methods."
        [arXiv: 1806.08693](https://arXiv.org/abs/1806.08693)

    Efficient implementation (and optimized controller) developed by

      - Hendrik Ranocha, Lisandro Dalcin, Matteo Parsani, David I. Ketcheson (2021)
        Optimized Runge-Kutta Methods with Automatic Step Size Control for
        Compressible Computational Fluid Dynamics
        [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)"""
)
Base.@kwdef struct SSPRK43{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, four-stage explicit strong stability preserving (SSP) method.",
    "SSPRK432",
    references = "Gottlieb, Sigal, David I. Ketcheson, and Chi-Wang Shu.
    Strong stability preserving Runge-Kutta and multistep time discretizations.
    World Scientific, 2011.
    Example 6.1"
)
Base.@kwdef struct SSPRK432{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A second-order, three-step explicit strong stability preserving (SSP) linear multistep method.
This method does not come with an error estimator and requires a fixed time step
size.",
    "SSPRKMSVS32",
    references = "Shu, Chi-Wang.
    Total-variation-diminishing time discretizations.
    SIAM Journal on Scientific and Statistical Computing 9, no. 6 (1988): 1073-1084.
    [DOI: 10.1137/0909073](https://doi.org/10.1137/0909073)"
)
Base.@kwdef struct SSPRKMSVS32{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A fourth-order, five-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK54",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207."
)
Base.@kwdef struct SSPRK54{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, five-stage explicit strong stability preserving (SSP) low-storage method.
Fixed timestep only.",
    "SSPRK53_2N1",
    references = "Higueras and T. Roldán.
    New third order low-storage SSP explicit Runge–Kutta methods
    arXiv:1809.04807v1."
)
Base.@kwdef struct SSPRK53_2N1{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A fourth-order, ten-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK104",
    references = "Ketcheson, David I.
    Highly efficient strong stability-preserving Runge–Kutta methods with
    low-storage implementations.
    SIAM Journal on Scientific Computing 30.4 (2008): 2113-2136."
)
Base.@kwdef struct SSPRK104{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, nine-stage explicit strong stability preserving (SSP) method.

Consider using `SSPRK43` instead, which uses the same main method and an
improved embedded method.",
    "SSPRK932",
    references = "Gottlieb, Sigal, David I. Ketcheson, and Chi-Wang Shu.
    Strong stability preserving Runge-Kutta and multistep time discretizations.
    World Scientific, 2011."
)
Base.@kwdef struct SSPRK932{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, four-step explicit strong stability preserving (SSP) linear multistep method.
This method does not come with an error estimator and requires a fixed time step
size.",
    "SSPRKMSVS43",
    references = "Shu, Chi-Wang.
    Total-variation-diminishing time discretizations.
    SIAM Journal on Scientific and Statistical Computing 9, no. 6 (1988): 1073-1084.
    [DOI: 10.1137/0909073](https://doi.org/10.1137/0909073)"
)
Base.@kwdef struct SSPRKMSVS43{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, seven-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK73",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207"
)
Base.@kwdef struct SSPRK73{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, five-stage explicit strong stability preserving (SSP) low-storage method.
Fixed timestep only.",
    "SSPRK53_H",
    references = "Higueras and T. Roldán.
    New third order low-storage SSP explicit Runge–Kutta methods
    arXiv:1809.04807v1."
)
Base.@kwdef struct SSPRK53_H{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "A third-order, three-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK33",
    references = "Shu, Chi-Wang, and Stanley Osher.
    Efficient implementation of essentially non-oscillatory shock-capturing schemes.
    Journal of Computational Physics 77.2 (1988): 439-471.
    https://doi.org/10.1016/0021-9991(88)90177-5"
)
Base.@kwdef struct SSPRK33{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "Optimal strong-stability-preserving Runge-Kutta time discretizations for discontinuous Galerkin methods",
    "KYKSSPRK42",
    references = "@article{kubatko2014optimal,
    title={Optimal strong-stability-preserving Runge--Kutta time discretizations for discontinuous Galerkin methods},
    author={Kubatko, Ethan J and Yeager, Benjamin A and Ketcheson, David I},
    journal={Journal of Scientific Computing},
    volume={60},
    pages={313--344},
    year={2014},
    publisher={Springer}}"
)
Base.@kwdef struct KYKSSPRK42{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "Optimal strong-stability-preserving Runge-Kutta time discretizations for discontinuous Galerkin methods",
    "KYK2014DGSSPRK_3S2",
    references = """@article{kubatko2014optimal,
    title={Optimal strong-stability-preserving Runge--Kutta time discretizations for discontinuous Galerkin methods},
    author={Kubatko, Ethan J and Yeager, Benjamin A and Ketcheson, David I},
    journal={Journal of Scientific Computing},
    volume={60},
    pages={313--344},
    year={2014},
    publisher={Springer}}"""
)
Base.@kwdef struct KYK2014DGSSPRK_3S2{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "Parametric relaxation Runge-Kutta method based on Heun's method (SSPRK(2,2)).
The parametric stabilization adds a term κ(u - u) to the ODE, modifying
Shu-Osher coefficients to enhance stability for large time steps.
SSP coefficient C = 1. The parameter `kappa` controls the stabilization strength;
when `kappa = 0`, this reduces to standard SSPRK22.
Fixed timestep only.",
    "pRRK22",
    references = """@article{liu2023high,
    title={High-order, large time-stepping integrators for scalar hyperbolic conservation laws},
    author={Liu, Lele and Burchard, Hans and Kuzmin, Dmitri and Shang, Sijun},
    year={2023},
    note={SSRN preprint 4401014}}"""
)
Base.@kwdef struct pRRK22{T, StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    kappa::T = 0.0
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "Parametric relaxation Runge-Kutta method based on RK(3,3) (SSPRK(3,3)).
The parametric stabilization adds a term κ(u - u) to the ODE, modifying
Shu-Osher coefficients to enhance stability for large time steps.
SSP coefficient C = 1. The parameter `kappa` controls the stabilization strength;
when `kappa = 0`, this reduces to standard SSPRK33.
Fixed timestep only.",
    "pRRK33",
    references = """@article{liu2023high,
    title={High-order, large time-stepping integrators for scalar hyperbolic conservation laws},
    author={Liu, Lele and Burchard, Hans and Kuzmin, Dmitri and Shang, Sijun},
    year={2023},
    note={SSRN preprint 4401014}}"""
)
Base.@kwdef struct pRRK33{T, StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    kappa::T = 0.0
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

@doc explicit_rk_docstring(
    "Parametric relaxation Runge-Kutta method based on Spiteri-Ruuth RK(5,4) (SSPRK(5,4)).
The parametric stabilization adds a term κ(u - u) to the ODE, modifying
Shu-Osher coefficients to enhance stability for large time steps.
SSP coefficient C ≈ 1.508. The parameter `kappa` controls the stabilization strength;
when `kappa = 0`, this reduces to standard SSPRK54.
Fixed timestep only.",
    "pRRK54",
    references = """@article{liu2023high,
    title={High-order, large time-stepping integrators for scalar hyperbolic conservation laws},
    author={Liu, Lele and Burchard, Hans and Kuzmin, Dmitri and Shang, Sijun},
    year={2023},
    note={SSRN preprint 4401014}}"""
)
Base.@kwdef struct pRRK54{T, StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    kappa::T = 0.0
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end
