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
- Di Marzo G. RODAS5(4) – Méthodes de Rosenbrock d'ordre 5(4) adaptées aux problèmes
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

=#

# Documentation for Rosenbrock methods with step_limiter

@doc rosenbrock_wolfbrandt_docstring(
"""
An Order 2/3 L-Stable Rosenbrock-W method which is good for very stiff equations with oscillations at low tolerances. 2nd order stiff-aware interpolation.
""",
"Rosenbrock23",
references = """
- Shampine L.F. and Reichelt M., (1997) The MATLAB ODE Suite, SIAM Journal of
  Scientific Computing, 18 (1), pp. 1-22.
""",
with_step_limiter = true) Rosenbrock23

@doc rosenbrock_wolfbrandt_docstring(
"""
An Order 3/2 A-Stable Rosenbrock-W method which is good for mildly stiff equations without oscillations at low tolerances. Note that this method is prone to instability in the presence of oscillations, so use with caution. 2nd order stiff-aware interpolation.
""",
"Rosenbrock32",
references = """
- Shampine L.F. and Reichelt M., (1997) The MATLAB ODE Suite, SIAM Journal of
  Scientific Computing, 18 (1), pp. 1-22.
""",
with_step_limiter = true) Rosenbrock32

@doc rosenbrock_docstring(
"""
3rd order A-stable and stiffly stable Rosenbrock method. Keeps high accuracy on discretizations of nonlinear parabolic PDEs.
""",
"ROS3P",
references = """
- Lang, J. & Verwer, ROS3P—An Accurate Third-Order Rosenbrock Solver Designed for
  Parabolic Problems J. BIT Numerical Mathematics (2001) 41: 731. doi:10.1023/A:1021900219772
""",
with_step_limiter = true) ROS3P

@doc rosenbrock_docstring(
"""
3rd order A-stable and stiffly stable Rosenbrock method.
""",
"Rodas3",
references = """
- Sandu, Verwer, Van Loon, Carmichael, Potra, Dabdub, Seinfeld, Benchmarking stiff ode solvers for atmospheric chemistry problems-I. 
  implicit vs explicit, Atmospheric Environment, 31(19), 3151-3166, 1997.
""",
with_step_limiter=true) Rodas3

@doc rosenbrock_wolfbrandt_docstring(
"""
An Order 2/3 L-Stable Rosenbrock-W method for stiff ODEs and DAEs in mass matrix form. 2nd order stiff-aware interpolation and additional error test for interpolation.
""",
"Rodas23W",
references = """
- Steinebach G., Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -
  Preprint 2024. Proceedings of the JuliaCon Conferences.
  https://proceedings.juliacon.org/papers/eb04326e1de8fa819a3595b376508a40
""",
with_step_limiter = true) Rodas23W

@doc rosenbrock_docstring(
"""
3rd order A-stable and stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant
and additional error test for interpolation. Keeps accuracy on discretizations of linear parabolic PDEs.
""",
"Rodas3P",
references = """
- Steinebach G., Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -
  Preprint 2024. Proceedings of the JuliaCon Conferences.
  https://proceedings.juliacon.org/papers/eb04326e1de8fa819a3595b376508a40
""",
with_step_limiter=true) Rodas3P

@doc rosenbrock_docstring(
"""
A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant
""",
"Rodas4",
references = """
- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
""",
with_step_limiter=true) Rodas4

@doc rosenbrock_docstring(
"""
A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant
""",
"Rodas42",
references = """
- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
""",
with_step_limiter=true) Rodas42

@doc rosenbrock_docstring(
"""
4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant. 4th order
on linear parabolic problems and 3rd order accurate on nonlinear parabolic problems (as opposed to
lower if not corrected).
""",
"Rodas4P",
references = """
- Steinebach, G., Rentrop, P., An adaptive method of lines approach for modelling flow and transport in rivers. 
  Adaptive method of lines , Wouver, A. Vande, Sauces, Ph., Schiesser, W.E. (ed.),S. 181-205,Chapman & Hall/CRC, 2001,
- Steinebach, G., Order-reduction of ROW-methods for DAEs and method of lines  applications. 
  Preprint-Nr. 1741, FB Mathematik, TH Darmstadt, 1995.
""",
with_step_limiter=true) Rodas4P

@doc rosenbrock_wolfbrandt_docstring(
"""
A 4th order L-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant. 4th order
on linear parabolic problems and 3rd order accurate on nonlinear parabolic problems. It is an improvement
of Rodas4P and in case of inexact Jacobians a second order W method.
""",
"Rodas4P2",
references = """
- Steinebach G., Improvement of Rosenbrock-Wanner Method RODASP, In: Reis T., Grundel S., Schöps S. (eds) 
  Progress in Differential-Algebraic Equations II. Differential-Algebraic Equations Forum. Springer, Cham., 165-184, 2020.
""",
with_step_limiter=true) Rodas4P2

@doc rosenbrock_docstring(
"""
A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware 4th order interpolant.
""",
"Rodas5",
references = """
- Di Marzo G. RODAS5(4) – Méthodes de Rosenbrock d'ordre 5(4) adaptées aux problèmes
  différentiels-algébriques. MSc mathematics thesis, Faculty of Science,
  University of Geneva, Switzerland.
""",
with_step_limiter=true) Rodas5

@doc rosenbrock_wolfbrandt_docstring(
"""
A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware 4th order interpolant.
Has improved stability in the adaptive time stepping embedding.
""",
"Rodas5P",
references = """
- Steinebach G. Construction of Rosenbrock–Wanner method Rodas5P and numerical benchmarks
  within the Julia Differential Equations package.
  In: BIT Numerical Mathematics, 63(2), 2023. doi:10.1007/s10543-023-00967-x
""",
with_step_limiter=true) Rodas5P

@doc rosenbrock_wolfbrandt_docstring(
"""
Variant of Rodas5P with additional residual control.
""",
"Rodas5Pr",
references = """
- Steinebach G. Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -
  Preprint 2024. Proceedings of the JuliaCon Conferences.
  https://proceedings.juliacon.org/papers/eb04326e1de8fa819a3595b376508a40
""",
with_step_limiter=true) Rodas5Pr

@doc rosenbrock_wolfbrandt_docstring(
"""
Variant of Rodas5P with modified embedded scheme.
""",
"Rodas5Pe",
references = """
- Steinebach G. Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -
  Preprint 2024. Proceedings of the JuliaCon Conferences.
  https://proceedings.juliacon.org/papers/eb04326e1de8fa819a3595b376508a40
""",
with_step_limiter=true) Rodas5Pe

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
    :Rodas5Pr
]
    @eval begin
        struct $Alg{CS, AD, F, P, FDT, ST, CJ, StepLimiter, StageLimiter} <:
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

"""
$(rosenbrock_wolfbrandt_docstring(
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

@doc rosenbrock_docstring(
"""
A 2nd order L-stable Rosenbrock method with 2 internal stages.
""",
"ROS2",
references = """
- J. G. Verwer et al. (1999): A second-order Rosenbrock method applied to photochemical dispersion problems
  https://doi.org/10.1137/S1064827597326651
""") ROS2

@doc rosenbrock_docstring(
"""
2nd order stiffly accurate Rosenbrock method with 3 internal stages with (Rinf=0).
For problems with medium stiffness the convergence behaviour is very poor and it is recommended to use
[`ROS2S`](@ref) instead.
""",
"ROS2PR",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") ROS2PR

@doc rosenbrock_wolfbrandt_docstring(
"""
2nd order stiffly accurate Rosenbrock-Wanner W-method with 3 internal stages with B_PR consistent of order 2 with (Rinf=0).
""",
"ROS2S",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") ROS2S

@doc rosenbrock_docstring(
"""
3rd order L-stable Rosenbrock method with 3 internal stages with an embedded strongly
A-stable 2nd order method.
""",
"ROS3",
references = """
- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
""") ROS3

@doc rosenbrock_docstring(
"""
3nd order stiffly accurate Rosenbrock method with 3 internal stages with B_PR consistent of order 3, which is strongly A-stable with Rinf~=-0.73.
""",
"ROS3PR",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") ROS3PR

@doc rosenbrock_docstring(
"""
3nd order stiffly accurate Rosenbrock method with 3 internal stages with B_PR consistent of order 3, which is strongly A-stable with Rinf~=-0.73.
Convergence with order 4 for the stiff case, but has a poor accuracy.
""",
"Scholz4_7",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") Scholz4_7

@doc rosenbrock_wolfbrandt_docstring(
"""
A 4th order L-stable Rosenbrock-W method.
""",
"ROS34PW1a",
references = """
- Rang, Joachim and Angermann, L (2005): New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1.
  BIT Numerical Mathematics, 45, 761--787.
""") ROS34PW1a

@doc rosenbrock_wolfbrandt_docstring(
"""
A 4th order L-stable Rosenbrock-W method.
""",
"ROS34PW1b",
references = """
- Rang, Joachim and Angermann, L (2005): New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1.
  BIT Numerical Mathematics, 45, 761--787.
""") ROS34PW1b

@doc rosenbrock_wolfbrandt_docstring(
"""
A 4th order stiffy accurate Rosenbrock-W method for PDAEs.
""",
"ROS34PW2",
references = """
- Rang, Joachim and Angermann, L (2005): New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1.
  BIT Numerical Mathematics, 45, 761--787.
""") ROS34PW2

@doc rosenbrock_wolfbrandt_docstring(
"""
A 4th order strongly A-stable (Rinf~0.63) Rosenbrock-W method.
""",
"ROS34PW3",
references = """
- Rang, Joachim and Angermann, L (2005): New Rosenbrock W-methods of order 3 for partial differential algebraic equations of index 1.
  BIT Numerical Mathematics, 45, 761--787.
""") ROS34PW3

@doc rosenbrock_wolfbrandt_docstring(
"""
3rd order stiffly accurate Rosenbrock-Wanner W-method with 4 internal stages,
B_PR consistent of order 2.
The order of convergence decreases if medium stiff problems are considered.
""",
"ROS34PRw",
references = """
- Joachim Rang, Improved traditional Rosenbrock–Wanner methods for stiff ODEs and DAEs,
  Journal of Computational and Applied Mathematics,
  https://doi.org/10.1016/j.cam.2015.03.010
""") ROS34PRw

@doc rosenbrock_docstring(
"""
3rd order stiffly accurate Rosenbrock method with 4 internal stages,
B_PR consistent of order 2 with Rinf=0.
The order of convergence decreases if medium stiff problems are considered, but it has good results for very stiff cases.
""",
"ROS3PRL",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") ROS3PRL

@doc rosenbrock_docstring(
"""
3rd order stiffly accurate Rosenbrock method with 4 internal stages,
B_PR consistent of order 3.
The order of convergence does NOT decreases if medium stiff problems are considered as it does for [`ROS3PRL`](@ref).
""",
"ROS3PRL2",
references = """
- Rang, Joachim (2014): The Prothero and Robinson example:
  Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
  https://doi.org/10.24355/dbbs.084-201408121139-0
""") ROS3PRL2

@doc rosenbrock_wolfbrandt_docstring(
"""
4rd order L-stable Rosenbrock-Krylov method with 4 internal stages,
with a 3rd order embedded method which is strongly A-stable with Rinf~=0.55. (when using exact Jacobians)
""",
"ROK4a",
references = """
- Tranquilli, Paul and Sandu, Adrian (2014):
  Rosenbrock--Krylov Methods for Large Systems of Differential Equations
  https://doi.org/10.1137/130923336
""") ROK4a

@doc rosenbrock_docstring(
"""
An A-stable 4th order Rosenbrock method.
""",
"RosShamp4",
references = """
- L. F. Shampine, Implementation of Rosenbrock Methods, ACM Transactions on
  Mathematical Software (TOMS), 8: 2, 93-113. doi:10.1145/355993.355994
""") RosShamp4

@doc rosenbrock_docstring(
"""
A 4th order D-stable Rosenbrock method.
""",
"Veldd4",
references = """
- van Veldhuizen, D-stability and Kaps-Rentrop-methods, M. Computing (1984) 32: 229.
  doi:10.1007/BF02243574
""") Veldd4

@doc rosenbrock_wolfbrandt_docstring(
"""
A 4th order A-stable Rosenbrock method.
""",
"Velds4",
references = """
- van Veldhuizen, D-stability and Kaps-Rentrop-methods, M. Computing (1984) 32: 229.
  doi:10.1007/BF02243574
""") Velds4

@doc rosenbrock_docstring(
"""
An efficient 4th order Rosenbrock method.
""",
"GRK4T",
references = """
- Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four with stepsize control
  for stiff ordinary differential equations. P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495
""") GRK4T

@doc rosenbrock_docstring(
"""
An A-stable 4th order Rosenbrock method. Essentially "anti-L-stable" but efficient.
""",
"GRK4A",
references = """
- Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four with stepsize control
  for stiff ordinary differential equations. P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495
""") GRK4A

@doc rosenbrock_docstring(
"""
A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware 3rd order interpolant
""",
"Ros4LStab",
references = """
- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)
""") Ros4LStab

for Alg in [
    :ROS2,
    :ROS2PR,
    :ROS2S,
    :ROS3,
    :ROS3PR,
    :Scholz4_7,
    :ROS34PW1a,
    :ROS34PW1b,
    :ROS34PW2,
    :ROS34PW3,
    :ROS34PRw,
    :ROS3PRL,
    :ROS3PRL2,
    :ROK4a,
    :RosShamp4,
    :Veldd4,
    :Velds4,
    :GRK4T,
    :GRK4A,
    :Ros4LStab
]
    @eval begin
        struct $Alg{CS, AD, F, P, FDT, ST, CJ} <:
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

