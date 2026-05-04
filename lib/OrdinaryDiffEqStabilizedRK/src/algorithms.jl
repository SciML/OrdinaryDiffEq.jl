@doc generic_solver_docstring(
    """High stability for real eigenvalues. Second order method. Exhibits high stability for real eigenvalues
    and is smoothened to allow for moderate sized complex eigenvalues.""",
    "ROCK2",
    "Stabilized Explicit Method.",
    """Assyr Abdulle, Alexei A. Medovikov. Second Order Chebyshev Methods based on Orthogonal Polynomials.
    Numerische Mathematik, 90 (1), pp 1-18, 2001. doi: https://dx.doi.org/10.1007/s002110100292""",
    """
    - `min_stages`: The minimum degree of the Chebyshev polynomial.
    - `max_stages`: The maximumdegree of the Chebyshev polynomial.
    - `eigen_est`: function of the form
        `(integrator) -> integrator.eigen_est = upper_bound`,
        where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
        If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
    """,
    """
    min_stages = 0,
    max_stages = 200,
    eigen_est = nothing,
    """
)
struct ROCK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    min_stages::Int
    max_stages::Int
    eigen_est::E
end
function ROCK2(; min_stages = 0, max_stages = 200, eigen_est = nothing)
    return ROCK2(min_stages, max_stages, eigen_est)
end

@doc generic_solver_docstring(
    """High stability for real eigenvalues. Fourth order method. Exhibits high stability for real eigenvalues
    and is smoothened to allow for moderate sized complex eigenvalues.""",
    "ROCK4",
    "Stabilized Explicit Method.",
    """Assyr Abdulle. Fourth Order Chebyshev Methods With Recurrence Relation. 2002 Society for
    Industrial and Applied Mathematics Journal on Scientific Computing, 23(6), pp 2041-2054, 2001.
    doi: https://doi.org/10.1137/S1064827500379549""",
    """
    - `min_stages`: The minimum degree of the Chebyshev polynomial.
    - `max_stages`: The maximumdegree of the Chebyshev polynomial.
    - `eigen_est`: function of the form
        `(integrator) -> integrator.eigen_est = upper_bound`,
        where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
        If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
    """,
    """
    min_stages = 0,
    max_stages = 152,
    eigen_est = nothing,
    """
)
struct ROCK4{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    min_stages::Int
    max_stages::Int
    eigen_est::E
end
function ROCK4(; min_stages = 0, max_stages = 152, eigen_est = nothing)
    return ROCK4(min_stages, max_stages, eigen_est)
end

# SERK methods

for Alg in [:ESERK4, :ESERK5, :RKC, :TSRKC3]
    @eval begin
        struct $Alg{E} <: OrdinaryDiffEqAdaptiveAlgorithm
            eigen_est::E
        end
        $Alg(; eigen_est = nothing) = $Alg(eigen_est)
    end
end

@doc generic_solver_docstring(
    """Second order method. Exhibits high stability for real eigenvalues.""",
    "RKC",
    "Stabilized Explicit Method.",
    """B. P. Sommeijer, L. F. Shampine, J. G. Verwer. RKC: An Explicit Solver for Parabolic PDEs,
    Journal of Computational and Applied Mathematics, 88(2), pp 315-326, 1998. doi:
    https://doi.org/10.1016/S0377-0427(97)00219-7""",
    """
    - `eigen_est`: function of the form
        `(integrator) -> integrator.eigen_est = upper_bound`,
        where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
        If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
    """,
    """
    eigen_est = nothing,
    """
)
function RKC end

@doc generic_solver_docstring(
    """Fourth order method. Exhibits high stability for real eigenvalues
    and is smoothened to allow for moderate sized complex eigenvalues.""",
    "ESERK4",
    "Stabilized Explicit Method.",
    """J. Martín-Vaquero, B. Kleefeld. Extrapolated stabilized explicit Runge-Kutta methods,
    Journal of Computational Physics, 326, pp 141-155, 2016. doi:
    https://doi.org/10.1016/j.jcp.2016.08.042.""",
    """
    - `eigen_est`: function of the form
        `(integrator) -> integrator.eigen_est = upper_bound`,
        where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
        If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
    """,
    """
    eigen_est = nothing,
    """
)
function ESERK4 end

@doc generic_solver_docstring(
    """Fifth order method. Exhibits high stability for real eigenvalues
    and is smoothened to allow for moderate sized complex eigenvalues.""",
    "ESERK5",
    "Stabilized Explicit Method.",
    """J. Martín-Vaquero, A. Kleefeld. ESERK5: A fifth-order extrapolated stabilized explicit Runge-Kutta method,
    Journal of Computational and Applied Mathematics, 356, pp 22-36, 2019. doi:
    https://doi.org/10.1016/j.cam.2019.01.040.""",
    """
    - `eigen_est`: function of the form
        `(integrator) -> integrator.eigen_est = upper_bound`,
        where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
        If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
    """,
    """
    eigen_est = nothing,
    """
)
function ESERK5 end

@doc generic_solver_docstring(
    """Second order method.""",
    "SERK2",
    "Stabilized Explicit Method.",
    """@article{kleefeld2013serk2v2,
    title={SERK2v2: A new second-order stabilized explicit Runge-Kutta method for stiff problems},
    author={Kleefeld, B and Martin-Vaquero, J},
    journal={Numerical Methods for Partial Differential Equations},
    volume={29},
    number={1},
    pages={170--185},
    year={2013},
    publisher={Wiley Online Library}}""",
    """
    - `controller`: TBD
    - `eigen_est`: function of the form
        `(integrator) -> integrator.eigen_est = upper_bound`,
        where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
        If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
    """,
    """
    eigen_est = nothing,
    """
)
struct SERK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    eigen_est::E
end
SERK2(; eigen_est = nothing) = SERK2(eigen_est)

@doc generic_solver_docstring(
    """Third order method. Exhibits high stability for real eigenvalues.""",
    "TSRKC3",
    "Two-step Stabilized Explicit Method.",
    """A. V. Moisa. Third order two-step Runge-Kutta-Chebyshev methods,
    Journal of Computational and Applied Mathematics, 457, pp 116291, 2025. doi:
    https://doi.org/10.1016/j.cam.2024.116291""",
    """
    - `eigen_est`: function of the form
        `(integrator) -> integrator.eigen_est = upper_bound`,
        where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
        If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
    """,
    """
    eigen_est = nothing,
    """
)
function TSRKC3 end

@doc generic_solver_docstring(
    """First-order super-time-stepping method based on shifted Legendre polynomials.
    Monotone and convex-monotone stable for parabolic operators. The stage count s
    is chosen adaptively from the spectral radius so that the superstep scales as s²
    times the explicit timestep.""",
    "RKL1",
    "Stabilized Explicit Method.",
    """C. D. Meyer, D. S. Balsara, T. D. Aslam. A stabilized Runge-Kutta-Legendre method
    for explicit super-time-stepping of parabolic and mixed equations.
    Journal of Computational Physics, 257, pp 594-626, 2014.
    doi: https://doi.org/10.1016/j.jcp.2013.08.021""",
    """
    - `min_stages`: Minimum number of stages s (must be odd, >= 3).
    - `max_stages`: Maximum number of stages s.
    - `eigen_est`: Optional function `(integrator) -> integrator.eigen_est = upper_bound`
        providing an upper bound on the spectral radius. If not provided, estimated
        by power iteration.
    """,
    """
    min_stages = 3,
    max_stages = 200,
    eigen_est = nothing,
    """
)
struct RKL1{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    min_stages::Int
    max_stages::Int
    eigen_est::E
end
function RKL1(; min_stages = 3, max_stages = 200, eigen_est = nothing)
    min_s = max(3, min_stages)
    min_s = isodd(min_s) ? min_s : min_s + 1
    max_s = isodd(max_stages) ? max_stages : max_stages - 1
    max_s = max(max_s, min_s)
    if max_s < min_s
        max_s = min_s
    end
    return RKL1(min_s, max_s, eigen_est)
end

@doc generic_solver_docstring(
    """Second-order super-time-stepping method based on shifted Legendre polynomials.
    Monotone and convex-monotone stable for parabolic operators with no spurious
    staircasing. The stage count s is chosen adaptively from the spectral radius
    so that the superstep scales as s² times the explicit timestep. Only odd values
    of s are used to ensure adequate damping of shortest wavelength modes.""",
    "RKL2",
    "Stabilized Explicit Method.",
    """C. D. Meyer, D. S. Balsara, T. D. Aslam. A stabilized Runge-Kutta-Legendre method
    for explicit super-time-stepping of parabolic and mixed equations.
    Journal of Computational Physics, 257, pp 594-626, 2014.
    doi: https://doi.org/10.1016/j.jcp.2013.08.021""",
    """
    - `min_stages`: Minimum number of stages s (must be odd, >= 3).
    - `max_stages`: Maximum number of stages s.
    - `eigen_est`: Optional function `(integrator) -> integrator.eigen_est = upper_bound`
        providing an upper bound on the spectral radius. If not provided, estimated
        by power iteration.
    """,
    """
    min_stages = 3,
    max_stages = 200,
    eigen_est = nothing,
    """
)
struct RKL2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    min_stages::Int
    max_stages::Int
    eigen_est::E
end
function RKL2(; min_stages = 3, max_stages = 200, eigen_est = nothing)
    min_s = max(3, min_stages)
    min_s = isodd(min_s) ? min_s : min_s + 1
    max_s = isodd(max_stages) ? max_stages : max_stages - 1
    max_s = max(max_s, min_s)
    return RKL2(min_s, max_s, eigen_est)
end
