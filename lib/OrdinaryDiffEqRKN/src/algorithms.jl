@doc generic_solver_docstring(
    "Method of order three, which minimizes the amount of evaluated functions in each step. Fixed time steps only.
Second order ODE should not depend on the first derivative.",
    "IRKN3",
    "Improved Runge-Kutta-Nyström method",
    "@article{rabiei2012numerical,
    title={Numerical Solution of Second-Order Ordinary Differential Equations by Improved Runge-Kutta Nystrom Method},
    author={Rabiei, Faranak and Ismail, Fudziah and Norazak, S and Emadi, Saeid},
    publisher={Citeseer}}", "", ""
)
struct IRKN3 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring(
    "4th order explicit Runge-Kutta-Nyström method. Allows acceleration to depend on velocity.",
    "Nystrom4",
    "Improved Runge-Kutta-Nyström method",
    "E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
    Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
    Springer-Verlag.", "", ""
)
struct Nystrom4 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring(
    "A 4th order explicit method which can be applied directly to second order ODEs.
In particular, this method allows the acceleration equation to depend on the velocity.",
    "FineRKN4",
    "Improved Runge-Kutta-Nyström method",
    "@article{fine1987low,
    title={Low order practical {R}unge-{K}utta-{N}ystr{\"o}m methods},
    author={Fine, Jerry Michael},
    journal={Computing},
    volume={38},
    number={4},
    pages={281--297},
    year={1987},
    publisher={Springer}}", "", ""
)
struct FineRKN4 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "A 5th order explicit method which can be applied directly to second order ODEs.
In particular, this method allows the acceleration equation to depend on the velocity.",
    "FineRKN5",
    "Improved Runge-Kutta-Nyström method",
    "@article{fine1987low,
    title={Low order practical {R}unge-{K}utta-{N}ystr{\"o}m methods},
    author={Fine, Jerry Michael},
    journal={Computing},
    volume={38},
    number={4},
    pages={281--297},
    year={1987},
    publisher={Springer}}", "", ""
)
struct FineRKN5 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "A 4th order explicit method.
Used directly on second order ODEs, where the acceleration is independent from velocity
(ODE Problem is not dependent on the first derivative).",
    "Nystrom4VelocityIndependent",
    "Improved Runge-Kutta-Nyström method",
    "E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
    Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
    Springer-Verlag.", "", ""
)
struct Nystrom4VelocityIndependent <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring(
    "Improves Runge-Kutta-Nyström method of order four,
    which minimizes the amount of evaluated functions in each step.
    Fixed time steps only.
    Second order ODE should not be dependent on the first derivative.
    Recommended for smooth problems with expensive functions to evaluate.",
    "IRKN4",
    "Improved Runge-Kutta-Nyström method",
    "@article{rabiei2012numerical,
    title={Numerical Solution of Second-Order Ordinary Differential Equations by Improved Runge-Kutta Nystrom Method},
    author={Rabiei, Faranak and Ismail, Fudziah and Norazak, S and Emadi, Saeid},
    publisher={Citeseer}}", "", ""
)
struct IRKN4 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring(
    "A 5th order explicit method.
Used directly on second order ODEs, where the acceleration is independent from velocity
(ODE Problem is not dependent on the first derivative).",
    "Nystrom5VelocityIndependent",
    "Improved Runge-Kutta-Nyström method",
    "E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
    Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
    Springer-Verlag.", "", ""
)
struct Nystrom5VelocityIndependent <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring(
    "4th order explicit method.
    The second order ODE should not depend on the first derivative.",
    "DPRKN4",
    "Improved Runge-Kutta-Nyström method",
    "@article{Dormand1987FamiliesOR,
    title={Families of Runge-Kutta-Nystrom Formulae},
    author={J. R. Dormand and Moawwad E. A. El-Mikkawy and P. J. Prince},
    journal={Ima Journal of Numerical Analysis},
    year={1987},
    volume={7},
    pages={235-250}}", "", ""
)
struct DPRKN4 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "5th order explicit method.
    The second order ODE should not depend on the first derivative.",
    "DPRKN5",
    "Improved Runge-Kutta-Nyström method",
    "@article{Bettis1973ARN,
    title={A Runge-Kutta Nystrom algorithm},
    author={Dale G. Bettis},
    journal={Celestial mechanics},
    year={1973},
    volume={8},
    pages={229-233},
    publisher={Springer}}", "", ""
)
struct DPRKN5 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "6th order explicit method.
The second order ODE should not depend on the first derivative. Free 6th order interpolant",
    "DPRKN6",
    "Improved Runge-Kutta-Nyström method",
    "@article{Dormand1987FamiliesOR,
    title={Families of Runge-Kutta-Nystrom Formulae},
    author={J. R. Dormand and Moawwad E. A. El-Mikkawy and P. J. Prince},
    journal={Ima Journal of Numerical Analysis},
    year={1987},
    volume={7},
    pages={235-250}}", "", ""
)
struct DPRKN6 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "6th order explicit method.
    The second order ODE should not depend on the first derivative.
    Compared to [`DPRKN6`](@ref), this method has smaller truncation error coefficients
    which leads to performance gain when only the main solution points are considered.",
    "DPRKN6FM",
    "Improved Runge-Kutta-Nyström method",
    "@article{Dormand1987FamiliesOR,
    title={Families of Runge-Kutta-Nystrom Formulae},
    author={J. R. Dormand and Moawwad E. A. El-Mikkawy and P. J. Prince},
    journal={Ima Journal of Numerical Analysis},
    year={1987},
    volume={7},
    pages={235-250}}", "", ""
)
struct DPRKN6FM <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "8th order explicit method.
    The second order ODE should not depend on the first derivative.
    Not as efficient as [`DPRKN12`](@ref) when high accuracy is needed,
    however this solver is competitive with [`DPRKN6`](@ref) at lax tolerances and,
    depending on the problem, might be a good option between performance and accuracy.",
    "DPRKN8",
    "Improved Runge-Kutta-Nyström method",
    "@article{dormand1987high,
    title={High-order embedded Runge-Kutta-Nystrom formulae},
    author={Dormand, JR and El-Mikkawy, MEA and Prince, PJ},
    journal={IMA Journal of Numerical Analysis},
    volume={7},
    number={4},
    pages={423--430},
    year={1987},
    publisher={Oxford University Press}}", "", ""
)
struct DPRKN8 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "12th order explicit method.
    The second order ODE should not depend on the first derivative.
    Most efficient when high accuracy is needed.",
    "DPRKN12",
    "Improved Runge-Kutta-Nyström method",
    "@article{dormand1987high,
    title={High-order embedded Runge-Kutta-Nystrom formulae},
    author={Dormand, JR and El-Mikkawy, MEA and Prince, PJ},
    journal={IMA Journal of Numerical Analysis},
    volume={7},
    number={4},
    pages={423--430},
    year={1987},
    publisher={Oxford University Press}}", "", ""
)
struct DPRKN12 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "Embedded 4(3) pair of explicit methods.
    Integrates the periodic properties of the harmonic oscillator exactly.
    The second order ODE should not depend on the first derivative.
    Uses adaptive step size control. This method is extra efficient on periodic problems.",
    "ERKN4",
    "Improved Runge-Kutta-Nyström method",
    "@article{demba2017embedded,
    title={An Embedded 4 (3) Pair of Explicit Trigonometrically-Fitted Runge-Kutta-Nystr{\"o}m Method for Solving Periodic Initial Value Problems},
    author={Demba, MA and Senu, N and Ismail, F},
    journal={Applied Mathematical Sciences},
    volume={11},
    number={17},
    pages={819--838},
    year={2017}}", "", ""
)
struct ERKN4 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "Embedded 5(4) pair of explicit methods.
    Integrates the periodic properties of the harmonic oscillator exactly.
    The second order ODE should not depend on the first derivative.
    Uses adaptive step size control. This method is extra efficient on periodic problems.",
    "ERKN5",
    "Improved Runge-Kutta-Nyström method",
    "@article{demba20165,
    title={A 5 (4) Embedded Pair of Explicit Trigonometrically-Fitted Runge--Kutta--Nystr{\"o}m Methods for the Numerical Solution of Oscillatory Initial Value     Problems},
    author={Demba, Musa A and Senu, Norazak and Ismail, Fudziah},
    journal={Mathematical and Computational Applications},
    volume={21},
    number={4},
    pages={46},
    year={2016},
    publisher={Multidisciplinary Digital Publishing Institute}}", "", ""
)
struct ERKN5 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "Embedded pair of explicit methods.
    Integrates the periodic properties of the harmonic oscillator exactly.
    The second order ODE should not depend on the first derivative.
    Uses adaptive step size control. This method is extra efficient on periodic problems.",
    "ERKN7",
    "Improved Runge-Kutta-Nyström method",
    "@article{SimosOnHO,
    title={On high order Runge-Kutta-Nystr{\"o}m pairs},
    author={Theodore E. Simos and Ch. Tsitouras},
    journal={J. Comput. Appl. Math.},
    volume={400},
    pages={113753}}", "", ""
)
struct ERKN7 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

@doc generic_solver_docstring(
    "3 stage fourth order method to solve second order linear inhomogeneous IVPs.
Does not include an adaptive method. Solves for for d-dimensional differential systems of second order linear inhomogeneous equations.

!!! warning
This method is only fourth order for these systems, the method is second order otherwise!",
    "RKN4",
    "Improved Runge-Kutta-Nyström method",
    "@article{MONTIJANO2024115533,
    title = {Explicit Runge–Kutta–Nyström methods for the numerical solution of second order linear inhomogeneous IVPs},
    author = {J.I. Montijano and L. Rández and M. Calvo},
    journal = {Journal of Computational and Applied Mathematics},
    volume = {438},
    pages = {115533},
    year = {2024},}", "", ""
)
struct RKN4 <: OrdinaryDiffEqAlgorithm end
