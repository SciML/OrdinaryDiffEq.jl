# Rosenbrock Methods

# for Rosenbrock methods with step_limiter
for (Alg, desc, refs, is_W) in [
        (:Rosenbrock23, "An Order 2/3 L-Stable Rosenbrock-W method which is good for very stiff equations with oscillations at low tolerances. 2nd order stiff-aware interpolation.", "- Shampine L.F. and Reichelt M., (1997) The MATLAB ODE Suite, SIAM Journal of\n  Scientific Computing, 18 (1), pp. 1-22.", true),
        (:Rosenbrock32, "An Order 3/2 A-Stable Rosenbrock-W method which is good for mildly stiff equations without oscillations at low tolerances. Note that this method is prone to instability in the presence of oscillations, so use with caution. 2nd order stiff-aware interpolation.", "- Shampine L.F. and Reichelt M., (1997) The MATLAB ODE Suite, SIAM Journal of\n  Scientific Computing, 18 (1), pp. 1-22.", true),
        (:ROS3P, "3rd order A-stable and stiffly stable Rosenbrock method. Keeps high accuracy on discretizations of nonlinear parabolic PDEs.", "- Lang, J. & Verwer, ROS3P—An Accurate Third-Order Rosenbrock Solver Designed for\n  Parabolic Problems J. BIT Numerical Mathematics (2001) 41: 731. doi:10.1023/A:1021900219772", false),
        (:Rodas3, "3rd order A-stable and stiffly stable Rosenbrock method.", "- Sandu, Verwer, Van Loon, Carmichael, Potra, Dabdub, Seinfeld, Benchmarking stiff ode solvers for atmospheric chemistry problems-I. \n  implicit vs explicit, Atmospheric Environment, 31(19), 3151-3166, 1997.", false),
        (:Rodas23W, "An Order 2/3 L-Stable Rosenbrock-W method for stiff ODEs and DAEs in mass matrix form. 2nd order stiff-aware interpolation and additional error test for interpolation.", "- Steinebach G., Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -\n  Preprint 2024. Proceedings of the JuliaCon Conferences.\n  https://proceedings.juliacon.org/papers/eb04326e1de8fa819a3595b376508a40", true),
        (:Rodas3P, "3rd order A-stable and stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant\nand additional error test for interpolation. Keeps accuracy on discretizations of linear parabolic PDEs.", "- Steinebach G., Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -\n  Preprint 2024. Proceedings of the JuliaCon Conferences.\n  https://proceedings.juliacon.org/papers/eb04326e1de8fa819a3595b376508a40", false),
        (:Rodas4, "A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant", "- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and\n  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)", false),
        (:Rodas42, "A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant", "- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and\n  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)", false),
        (:Rodas4P, "4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant. 4th order\non linear parabolic problems and 3rd order accurate on nonlinear parabolic problems (as opposed to\nlower if not corrected).", "- Steinebach, G., Rentrop, P., An adaptive method of lines approach for modelling flow and transport in rivers. \n  Adaptive method of lines , Wouver, A. Vande, Sauces, Ph., Schiesser, W.E. (ed.),S. 181-205,Chapman & Hall/CRC, 2001,\n- Steinebach, G., Order-reduction of ROW-methods for DAEs and method of lines  applications. \n  Preprint-Nr. 1741, FB Mathematik, TH Darmstadt, 1995.", false),
        (:Rodas4P2, "A 4th order L-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant. 4th order\non linear parabolic problems and 3rd order accurate on nonlinear parabolic problems. It is an improvement\nof Rodas4P and in case of inexact Jacobians a second order W method.", "- Steinebach G., Improvement of Rosenbrock-Wanner Method RODASP, In: Reis T., Grundel S., Schöps S. (eds) \n  Progress in Differential-Algebraic Equations II. Differential-Algebraic Equations Forum. Springer, Cham., 165-184, 2020.", true),
        (:Rodas5, "A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware 4th order interpolant.", "- Di Marzo G. RODAS5(4) – Méthodes de Rosenbrock d'ordre 5(4) adaptées aux problèmes\n  différentiels-algébriques. MSc mathematics thesis, Faculty of Science,\n  University of Geneva, Switzerland.", false),
        (:Rodas5P, "A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware 4th order interpolant.\nHas improved stability in the adaptive time stepping embedding.", "- Steinebach G. Construction of Rosenbrock–Wanner method Rodas5P and numerical benchmarks\n  within the Julia Differential Equations package.\n  In: BIT Numerical Mathematics, 63(2), 2023. doi:10.1007/s10543-023-00967-x", true),
        (:Rodas5Pe, "Variant of Rodas5P with modified embedded scheme.", "- Steinebach G. Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -\n  Preprint 2024. Proceedings of the JuliaCon Conferences.\n  https://proceedings.juliacon.org/papers/eb04326e1de8fa819a3595b376508a40", true),
        (:Rodas5Pr, "Variant of Rodas5P with additional residual control.", "- Steinebach G. Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -\n  Preprint 2024. Proceedings of the JuliaCon Conferences.\n  https://proceedings.juliacon.org/papers/eb04326e1de8fa819a3595b376508a40", true),
        (:Rodas6P, "A 6th order A-stable stiffly stable Rosenbrock method with a stiff-aware 5th order interpolant.", "- Steinebach G. Construction of Rosenbrock–Wanner method Rodas6P.\n  to prepare, 2025", true),
    ]
    @eval begin
        @doc $(is_W ? rosenbrock_wolfbrandt_docstring(desc, String(Alg), references = refs, with_step_limiter = true) : rosenbrock_docstring(desc, String(Alg), references = refs, with_step_limiter = true)) struct $Alg{CS, AD, F, P, FDT, ST, CJ, StepLimiter, StageLimiter} <:
            OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
            linsolve::F
            precs::P
            step_limiter!::StepLimiter
            stage_limiter!::StageLimiter
            autodiff::AD
        end
        function $Alg(;
                chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
                standardtag = Val{true}(), concrete_jac = nothing,
                diff_type = Val{:forward}(), linsolve = nothing,
                precs = DEFAULT_PRECS, step_limiter! = trivial_limiter!,
                stage_limiter! = trivial_limiter!
            )
            AD_choice, chunk_size, diff_type = _process_AD_choice(
                autodiff, chunk_size, diff_type
            )
            return $Alg{
                _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
                typeof(precs), diff_type, _unwrap_val(standardtag),
                _unwrap_val(concrete_jac), typeof(step_limiter!),
                typeof(stage_limiter!),
            }(
                linsolve, precs, step_limiter!,
                stage_limiter!, AD_choice
            )
        end
    end
end

"""
$(
    rosenbrock_wolfbrandt_docstring(
        """
        A 4th order L-stable Rosenbrock-W method (fixed step only).
        """,
        "RosenbrockW6S4OS",
        references = """
        https://doi.org/10.1016/j.cam.2009.09.017
        """
    )
)
"""
struct RosenbrockW6S4OS{CS, AD, F, P, FDT, ST, CJ} <:
    OrdinaryDiffEqRosenbrockAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    autodiff::AD
end
function RosenbrockW6S4OS(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing,
        precs = DEFAULT_PRECS
    )
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    return RosenbrockW6S4OS{
        _unwrap_val(chunk_size),
        typeof(AD_choice), typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag), _unwrap_val(concrete_jac),
    }(
        linsolve,
        precs, AD_choice
    )
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
        (:Ros4LStab, "A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant", "- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and\n  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)", false),
    ]
    @eval begin
        @doc $(is_W ? rosenbrock_wolfbrandt_docstring(desc, String(Alg), references = refs, with_step_limiter = false) : rosenbrock_docstring(desc, String(Alg), references = refs, with_step_limiter = false)) struct $Alg{CS, AD, F, P, FDT, ST, CJ} <:
            OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
            linsolve::F
            precs::P
            autodiff::AD
        end
        function $Alg(;
                chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
                standardtag = Val{true}(), concrete_jac = nothing,
                diff_type = Val{:forward}(), linsolve = nothing, precs = DEFAULT_PRECS
            )
            AD_choice, chunk_size, diff_type = _process_AD_choice(
                autodiff, chunk_size, diff_type
            )

            return $Alg{
                _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
                typeof(precs), diff_type, _unwrap_val(standardtag),
                _unwrap_val(concrete_jac),
            }(
                linsolve,
                precs, AD_choice
            )
        end
    end
end

################################################################################
# HybridExplicitImplicitRK — generic tableau-based hybrid explicit/linear-implicit method
################################################################################

struct HybridExplicitImplicitRK{TabType, CS, AD, F, P, FDT, ST, CJ, StepLimiter, StageLimiter} <:
    OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    tab::TabType
    order::Int
    linsolve::F
    precs::P
    step_limiter!::StepLimiter
    stage_limiter!::StageLimiter
    autodiff::AD
end

function HybridExplicitImplicitRK(
        tab;
        order,
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(), linsolve = nothing,
        precs = DEFAULT_PRECS, step_limiter! = trivial_limiter!,
        stage_limiter! = trivial_limiter!
    )
    AD_choice, chunk_size, diff_type = _process_AD_choice(
        autodiff, chunk_size, diff_type
    )
    return HybridExplicitImplicitRK{
        typeof(tab), _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!),
        typeof(stage_limiter!),
    }(
        tab, order, linsolve, precs, step_limiter!,
        stage_limiter!, AD_choice
    )
end

# Keyword-only constructor for remake support
function HybridExplicitImplicitRK(;
        tab,
        order,
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(), linsolve = nothing,
        precs = DEFAULT_PRECS, step_limiter! = trivial_limiter!,
        stage_limiter! = trivial_limiter!
    )
    return HybridExplicitImplicitRK(
        tab;
        order, chunk_size, autodiff, standardtag, concrete_jac,
        diff_type, linsolve, precs, step_limiter!, stage_limiter!
    )
end

"""
A 12-stage order 5(4) hybrid explicit/linear-implicit method for semi-explicit index-1 DAEs.
Differential variables are treated explicitly (like Tsit5), algebraic variables use Rosenbrock-type
linear-implicit stages. Only the small algebraic Jacobian block needs factorization.
For pure ODEs (no algebraic constraints), reduces to an explicit Runge-Kutta method.

References:
- Steinebach G., Rodas6P and Tsit5DA - two new Rosenbrock-type methods for DAEs.
  arXiv:2511.21252, 2025.
"""
Tsit5DA(; kwargs...) = HybridExplicitImplicitRK(Tsit5DATableau; order = 5, kwargs...)
