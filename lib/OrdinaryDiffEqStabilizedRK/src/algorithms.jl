
"""
Assyr Abdulle, Alexei A. Medovikov. Second Order Chebyshev Methods based on Orthogonal Polynomials.
Numerische Mathematik, 90 (1), pp 1-18, 2001. doi: https://dx.doi.org/10.1007/s002110100292

ROCK2: Stabilized Explicit Method.
Second order stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.

This method takes optional keyword arguments `min_stages`, `max_stages`, and `eigen_est`.
The function `eigen_est` should be of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix. If `eigen_est`
is not provided, `upper_bound` will be estimated using the power iteration.
"""
struct ROCK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    min_stages::Int
    max_stages::Int
    eigen_est::E
end
function ROCK2(; min_stages = 0, max_stages = 200, eigen_est = nothing)
    ROCK2(min_stages, max_stages, eigen_est)
end

"""
    ROCK4(; min_stages = 0, max_stages = 152, eigen_est = nothing)

Assyr Abdulle. Fourth Order Chebyshev Methods With Recurrence Relation. 2002 Society for
Industrial and Applied Mathematics Journal on Scientific Computing, 23(6), pp 2041-2054, 2001.
doi: https://doi.org/10.1137/S1064827500379549

ROCK4: Stabilized Explicit Method.
Fourth order stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.

This method takes optional keyword arguments `min_stages`, `max_stages`, and `eigen_est`.
The function `eigen_est` should be of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix. If `eigen_est`
is not provided, `upper_bound` will be estimated using the power iteration.
"""
struct ROCK4{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    min_stages::Int
    max_stages::Int
    eigen_est::E
end
function ROCK4(; min_stages = 0, max_stages = 152, eigen_est = nothing)
    ROCK4(min_stages, max_stages, eigen_est)
end

# SERK methods

for Alg in [:ESERK4, :ESERK5, :RKC]
    @eval begin
        struct $Alg{E} <: OrdinaryDiffEqAdaptiveAlgorithm
            eigen_est::E
        end
        $Alg(; eigen_est = nothing) = $Alg(eigen_est)
    end
end

"""
    RKC(; eigen_est = nothing)

B. P. Sommeijer, L. F. Shampine, J. G. Verwer. RKC: An Explicit Solver for Parabolic PDEs,
Journal of Computational and Applied Mathematics, 88(2), pp 315-326, 1998. doi:
https://doi.org/10.1016/S0377-0427(97)00219-7

RKC: Stabilized Explicit Method.
Second order stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues.

This method takes the keyword argument `eigen_est` of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix. If `eigen_est`
is not provided, `upper_bound` will be estimated using the power iteration.
"""
function RKC end

"""
    ESERK4(; eigen_est = nothing)

J. Martín-Vaquero, B. Kleefeld. Extrapolated stabilized explicit Runge-Kutta methods,
Journal of Computational Physics, 326, pp 141-155, 2016. doi:
https://doi.org/10.1016/j.jcp.2016.08.042.

ESERK4: Stabilized Explicit Method.
Fourth order extrapolated stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.

This method takes the keyword argument `eigen_est` of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
"""
function ESERK4 end

"""
    ESERK5(; eigen_est = nothing)

J. Martín-Vaquero, A. Kleefeld. ESERK5: A fifth-order extrapolated stabilized explicit Runge-Kutta method,
Journal of Computational and Applied Mathematics, 356, pp 22-36, 2019. doi:
https://doi.org/10.1016/j.cam.2019.01.040.

ESERK5: Stabilized Explicit Method.
Fifth order extrapolated stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.

This method takes the keyword argument `eigen_est` of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
"""
function ESERK5 end

struct SERK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    controller::Symbol
    eigen_est::E
end
SERK2(; controller = :PI, eigen_est = nothing) = SERK2(controller, eigen_est)

struct IRKC{CS, AD, F, F2, P, FDT, ST, CJ, K, T, E} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    controller::Symbol
    eigen_est::E
end

function IRKC(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, controller = :Standard, eigen_est = nothing)
    IRKC{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(eigen_est)}(linsolve, nlsolve, precs, κ, tol,
        extrapolant, controller, eigen_est)
end