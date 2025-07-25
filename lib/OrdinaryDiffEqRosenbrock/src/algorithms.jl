# Rosenbrock Methods

#=
#### Rosenbrock23, Rosenbrock32, ode23s, ModifiedRosenbrockIntegrator

- Shampine L.F. and Reichelt M., (1997) The MATLAB ODE Suite, SIAM Journal of
Scientific Computing, 18 (1), pp. 1-22.

#### ROS2

- J. G. Verwer et al. (1999): A second-order Rosenbrock method applied to photochemical dispersion problems
  https://doi.org/10.1137/S1064827597326651

#### ROS3P

- Lang, J. & Verwer, ROS3P—An Accurate Third-Order Rosenbrock Solver Designed for
  Parabolic Problems J. BIT Numerical Mathematics (2001) 41: 731. doi:10.1023/A:1021900219772

#### ROS3, Rodas3, Ros4LStab, Rodas4, Rodas42

- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)

#### ROS2PR, ROS2S, ROS3PR, Scholz4_7
-Rang, Joachim (2014): The Prothero and Robinson example:
 Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
 https://doi.org/10.24355/dbbs.084-201408121139-0

#### RosShamp4

- L. F. Shampine, Implementation of Rosenbrock Methods, ACM Transactions on
  Mathematical Software (TOMS), 8: 2, 93-113. doi:10.1145/355993.355994

#### Veldd4, Velds4

- van Veldhuizen, D-stability and Kaps-Rentrop-methods, M. Computing (1984) 32: 229.
  doi:10.1007/BF02243574

#### GRK4T, GRK4A

- Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four with stepsize control
  for stiff ordinary differential equations. P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495

#### Rodas23W, Rodas3P

- Steinebach G., Rodas23W / Rodas32P - a Rosenbrock-type method for DAEs with additional error estimate for dense output and Julia implementation,
 in progress

#### Rodas4P

- Steinebach G. Order-reduction of ROW-methods for DAEs and method of lines
  applications. Preprint-Nr. 1741, FB Mathematik, TH Darmstadt; 1995.

#### Rodas4P2
- Steinebach G. (2020) Improvement of Rosenbrock-Wanner Method RODASP.
  In: Reis T., Grundel S., Schoeps S. (eds) Progress in Differential-Algebraic Equations II.
  Differential-Algebraic Equations Forum. Springer, Cham. https://doi.org/10.1007/978-3-030-53905-4_6

#### Rodas5
- Di Marzo G. RODAS5(4) – Méthodes de Rosenbrock d’ordre 5(4) adaptées aux problemes
différentiels-algébriques. MSc mathematics thesis, Faculty of Science,
University of Geneva, Switzerland.

#### ROS34PRw
-Joachim Rang, Improved traditional Rosenbrock–Wanner methods for stiff ODEs and DAEs,
 Journal of Computational and Applied Mathematics,
 https://doi.org/10.1016/j.cam.2015.03.010

#### ROS3PRL, ROS3PRL2
-Rang, Joachim (2014): The Prothero and Robinson example:
 Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
 https://doi.org/10.24355/dbbs.084-201408121139-0

#### ROK4a
- Tranquilli, Paul and Sandu, Adrian (2014):
  Rosenbrock--Krylov Methods for Large Systems of Differential Equations
  https://doi.org/10.1137/130923336

#### Rodas5P
- Steinebach G.   Construction of Rosenbrock–Wanner method Rodas5P and numerical benchmarks within the Julia Differential Equations package.
 In: BIT Numerical Mathematics, 63(2), 2023

 #### Rodas23W, Rodas3P, Rodas5Pe, Rodas5Pr
- Steinebach G. Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -
 Preprint 2024
 https://github.com/hbrs-cse/RosenbrockMethods/blob/main/paper/JuliaPaper.pdf

 #### Rodas6P
- Steinebach G.   Construction of Rosenbrock–Wanner method Rodas6P , to prepare

=#

# for Rosenbrock methods with step_limiter
for Alg in [
    :Rosenbrock23,
    :Rosenbrock32,
    :ROS3P,
    :Rodas3,
    :Rodas23W,
    :Rodas3P,
    :Rodas4,
    :Rodas42,
    :Rodas4P,
    :Rodas4P2,
    :Rodas5,
    :Rodas5P,
    :Rodas5Pe,
    :Rodas5Pr,
    :Rodas6P]
    @eval begin
        @doc $(is_W ? rosenbrock_wolfbrandt_docstring(desc, String(Alg), references = refs, with_step_limiter = true) : rosenbrock_docstring(desc, String(Alg), references = refs, with_step_limiter = true)) struct $Alg{CS, AD, F, P, FDT, ST, CJ, StepLimiter, StageLimiter} <:
               OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
            linsolve::F
            precs::P
            step_limiter!::StepLimiter
            stage_limiter!::StageLimiter
            autodiff::AD
        end
        function $Alg(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
                standardtag = Val{true}(), concrete_jac = nothing,
                diff_type = Val{:forward}(), linsolve = nothing,
                precs = DEFAULT_PRECS, step_limiter! = trivial_limiter!,
                stage_limiter! = trivial_limiter!)
            AD_choice, chunk_size, diff_type = _process_AD_choice(
                autodiff, chunk_size, diff_type)
            $Alg{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
                typeof(precs), diff_type, _unwrap_val(standardtag),
                _unwrap_val(concrete_jac), typeof(step_limiter!),
                typeof(stage_limiter!)}(linsolve, precs, step_limiter!,
                stage_limiter!, AD_choice)
        end
    end
end

@doc rosenbrock_docstring(
    "An Order 2/3 L-Stable Rosenbrock-W method which is good for very stiff equations with oscillations at low tolerances. 2nd order stiff-aware interpolation. Strong stability for highly stiff systems. Good at high tolerances (>1e-2) for stiff problems. Recommended for highly stiff problems, systems with significant oscillations, low tolerance requirements.",
    "Rosenbrock23", with_step_limiter = true)
Rosenbrock23

@doc rosenbrock_docstring(
    "Efficient for medium tolerance stiff problems. A 5th order A-stable and stiffly stable embedded Rosenbrock method for differential-algebraic problems.",
    "Rodas5P", with_step_limiter = true)
Rodas5P

@doc rosenbrock_docstring(
    "Efficient for medium and high tolerance stiff problems. A 6th order A-stable and stiffly stable embedded Rosenbrock method for differential-algebraic problems.",
    "Rodas6P", with_step_limiter = true)
Rodas6P

@doc rosenbrock_docstring(
    "A 3/2-order L-stable Rosenbrock-W method optimized for stiff problems. Good balance of accuracy and computational efficiency.",
    "Rosenbrock32", with_step_limiter = true)
Rosenbrock32

@doc rosenbrock_docstring(
    "A 3rd-order accurate L-stable Rosenbrock method designed for parabolic problems. Particularly effective for reaction-diffusion equations.",
    "ROS3P", with_step_limiter = true)
ROS3P

@doc rosenbrock_docstring(
    "A 3rd-order accurate L-stable Rosenbrock method from Hairer and Wanner. Good general-purpose stiff ODE solver with moderate computational cost.",
    "Rodas3", with_step_limiter = true)
Rodas3

@doc rosenbrock_docstring(
    "A 4th-order accurate L-stable Rosenbrock method. Well-suited for moderately stiff problems with good efficiency.",
    "Rodas4", with_step_limiter = true)
Rodas4

@doc rosenbrock_docstring(
    "A 4th-order accurate L-stable Rosenbrock method with improved error estimation. Enhanced version of Rodas4 for better step size control.",
    "Rodas42", with_step_limiter = true)
Rodas42

@doc rosenbrock_docstring(
    "A 4th-order accurate L-stable Rosenbrock method designed for differential-algebraic equations (DAEs). Optimized for index-1 DAE problems.",
    "Rodas4P", with_step_limiter = true)
Rodas4P

@doc rosenbrock_docstring(
    "An improved 4th-order accurate L-stable Rosenbrock method for DAEs with enhanced stability properties.",
    "Rodas4P2", with_step_limiter = true)
Rodas4P2

@doc rosenbrock_docstring(
    "A 5th-order accurate L-stable Rosenbrock method for differential-algebraic problems. Higher accuracy but increased computational cost.",
    "Rodas5", with_step_limiter = true)
Rodas5

struct GeneralRosenbrock{CS, AD, F, ST, CJ, TabType} <:
       OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, Val{:forward}, ST, CJ}
    tableau::TabType
    factorization::F
    autodiff::AD
end

function GeneralRosenbrock(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        factorization = lu!, tableau = ROSENBROCK_DEFAULT_TABLEAU)
    AD_choice, chunk_size, diff_type = _process_AD_choice(
        autodiff, chunk_size, Val{:forward}())

    GeneralRosenbrock{
        _unwrap_val(chunk_size), typeof(AD_choice), typeof(factorization),
        _unwrap_val(standardtag), _unwrap_val(concrete_jac), typeof(tableau)}(tableau,
        factorization, AD_choice)
end

@doc rosenbrock_wolfbrandt_docstring(
    """
    A 4th order L-stable Rosenbrock-W method (fixed step only).
    """,
    "RosenbrockW6S4OS",
    references = """
    https://doi.org/10.1016/j.cam.2009.09.017
    """))
"""
struct RosenbrockW6S4OS{CS, AD, F, P, FDT, ST, CJ} <:
       OrdinaryDiffEqRosenbrockAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    autodiff::AD
end
function RosenbrockW6S4OS(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing,
        precs = DEFAULT_PRECS)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    RosenbrockW6S4OS{_unwrap_val(chunk_size),
        typeof(AD_choice), typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(linsolve,
        precs, AD_choice)
end

# Documentation for Rosenbrock methods without step_limiter

for (Alg, desc, refs, is_W) in [
    (:ROS2, "A 2nd order L-stable Rosenbrock method with 2 internal stages.", "- J. G. Verwer et al. (1999): A second-order Rosenbrock method applied to photochemical dispersion problems\n  https://doi.org/10.1137/S1064827597326651", false),
    (:ROS2PR, "2nd order stiffly accurate Rosenbrock method with 3 internal stages with (Rinf=0).\nFor problems with medium stiffness the convergence behaviour is very poor and it is recommended to use\n[`ROS2S`](@ref) instead.", "- Rang, Joachim (2014): The Prothero and Robinson example:\n  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.\n  https://doi.org/10.24355/dbbs.084-201408121139-0", false),
    (:ROS2S, "2nd order stiffly accurate Rosenbrock-Wanner W-method with 3 internal stages with B_PR consistent of order 2 with (Rinf=0).", "- Rang, Joachim (2014): The Prothero and Robinson example:\n  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.\n  https://doi.org/10.24355/dbbs.084-201408121139-0", true),
    (:ROS3, "3rd order L-stable Rosenbrock method with 3 internal stages with an embedded strongly\nA-stable 2nd order method.", "- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and\n  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)", false),
    (:ROS3PR, "3nd order stiffly accurate Rosenbrock method with 3 internal stages with B_PR consistent of order 3, which is strongly A-stable with Rinf~=-0.73.", "- Rang, Joachim (2014): The Prothero and Robinson example:\n  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.\n  https://doi.org/10.24355/dbbs.084-201408121139-0", false),
    (:Scholz4_7, "3nd order stiffly accurate Rosenbrock method with 3 internal stages with B_PR consistent of order 3, which is strongly A-stable with Rinf~=-0.73.\nConvergence with order 4 for the stiff case, but has a poor accuracy.", "- Rang, Joachim (2014): The Prothero and Robinson example:\n  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.\n  https://doi.org/10.24355/dbbs.084-201408121139-0", false),
    (:ROS34PW1a, "A 4th order L-stable Rosenbrock-W method.", "- Rang, Joachim and Angermann, L (2005): New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1.\n  BIT Numerical Mathematics, 45, 761--787.", true),
    (:ROS34PW1b, "A 4th order L-stable Rosenbrock-W method.", "- Rang, Joachim and Angermann, L (2005): New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1.\n  BIT Numerical Mathematics, 45, 761--787.", true),
    (:ROS34PW2, "A 4th order stiffy accurate Rosenbrock-W method for PDAEs.", "- Rang, Joachim and Angermann, L (2005): New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1.\n  BIT Numerical Mathematics, 45, 761--787.", true),
    (:ROS34PW3, "A 4th order strongly A-stable (Rinf~0.63) Rosenbrock-W method.", "- Rang, Joachim and Angermann, L (2005): New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1.\n  BIT Numerical Mathematics, 45, 761--787.", true),
    (:ROS34PRw, "3rd order stiffly accurate Rosenbrock-Wanner W-method with 4 internal stages,\nB_PR consistent of order 2.\nThe order of convergence decreases if medium stiff problems are considered.", "- Joachim Rang, Improved traditional Rosenbrock–Wanner methods for stiff ODEs and DAEs,\n  Journal of Computational and Applied Mathematics,\n  https://doi.org/10.1016/j.cam.2015.03.010", true),
    (:ROS3PRL, "3rd order stiffly accurate Rosenbrock method with 4 internal stages,\nB_PR consistent of order 2 with Rinf=0.\nThe order of convergence decreases if medium stiff problems are considered, but it has good results for very stiff cases.", "- Rang, Joachim (2014): The Prothero and Robinson example:\n  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.\n  https://doi.org/10.24355/dbbs.084-201408121139-0", false),
    (:ROS3PRL2, "3rd order stiffly accurate Rosenbrock method with 4 internal stages,\nB_PR consistent of order 3.\nThe order of convergence does NOT decreases if medium stiff problems are considered as it does for [`ROS3PRL`](@ref).", "- Rang, Joachim (2014): The Prothero and Robinson example:\n  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.\n  https://doi.org/10.24355/dbbs.084-201408121139-0", false),
    (:ROK4a, "4rd order L-stable Rosenbrock-Krylov method with 4 internal stages,\nwith a 3rd order embedded method which is strongly A-stable with Rinf~=0.55. (when using exact Jacobians)", "- Tranquilli, Paul and Sandu, Adrian (2014):\n  Rosenbrock--Krylov Methods for Large Systems of Differential Equations\n  https://doi.org/10.1137/130923336", true),
    (:RosShamp4, "An A-stable 4th order Rosenbrock method.", "- L. F. Shampine, Implementation of Rosenbrock Methods, ACM Transactions on\n  Mathematical Software (TOMS), 8: 2, 93-113. doi:10.1145/355993.355994", false),
    (:Veldd4, "A 4th order D-stable Rosenbrock method.", "- van Veldhuizen, D-stability and Kaps-Rentrop-methods, M. Computing (1984) 32: 229.\n  doi:10.1007/BF02243574", false),
    (:Velds4, "A 4th order A-stable Rosenbrock method.", "- van Veldhuizen, D-stability and Kaps-Rentrop-methods, M. Computing (1984) 32: 229.\n  doi:10.1007/BF02243574", true),
    (:GRK4T, "An efficient 4th order Rosenbrock method.", "- Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four with stepsize control\n  for stiff ordinary differential equations. P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495", false),
    (:GRK4A, "An A-stable 4th order Rosenbrock method. Essentially \"anti-L-stable\" but efficient.", "- Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four with stepsize control\n  for stiff ordinary differential equations. P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495", false),
    (:Ros4LStab, "A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant", "- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and\n  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)", false)
]
    @eval begin
        @doc $(is_W ? rosenbrock_wolfbrandt_docstring(desc, String(Alg), references = refs, with_step_limiter = false) : rosenbrock_docstring(desc, String(Alg), references = refs, with_step_limiter = false)) struct $Alg{CS, AD, F, P, FDT, ST, CJ} <:
               OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
            linsolve::F
            precs::P
            autodiff::AD
        end
        function $Alg(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
                standardtag = Val{true}(), concrete_jac = nothing,
                diff_type = Val{:forward}(), linsolve = nothing, precs = DEFAULT_PRECS)
            AD_choice, chunk_size, diff_type = _process_AD_choice(
                autodiff, chunk_size, diff_type)

            $Alg{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
                typeof(precs), diff_type, _unwrap_val(standardtag),
                _unwrap_val(concrete_jac)}(linsolve,
                precs, AD_choice)
        end
    end
end

