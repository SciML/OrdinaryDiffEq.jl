"""
    IRKN3

Improved Runge-Kutta-Nyström method of order three, which minimizes the amount of evaluated functions in each step. Fixed time steps only.

Second order ODE should not depend on the first derivative.

## References

@article{rabiei2012numerical,
title={Numerical Solution of Second-Order Ordinary Differential Equations by Improved Runge-Kutta Nystrom Method},
author={Rabiei, Faranak and Ismail, Fudziah and Norazak, S and Emadi, Saeid},
publisher={Citeseer}
}
"""
struct IRKN3 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    Nystrom4

A 4th order explicit Runge-Kutta-Nyström method which can be applied directly on second order ODEs. Can only be used with fixed time steps.

In case the ODE Problem is not dependent on the first derivative consider using
[`Nystrom4VelocityIndependent`](@ref) to increase performance.

## References

E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
Springer-Verlag.
"""
struct Nystrom4 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    FineRKN4()

A 4th order explicit Runge-Kutta-Nyström method which can be applied directly to second order ODEs.
In particular, this method allows the acceleration equation to depend on the velocity.

## References

```
@article{fine1987low,
  title={Low order practical {R}unge-{K}utta-{N}ystr{\"o}m methods},
  author={Fine, Jerry Michael},
  journal={Computing},
  volume={38},
  number={4},
  pages={281--297},
  year={1987},
  publisher={Springer}
}
```
"""
struct FineRKN4 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    FineRKN5()

A 5th order explicit Runge-Kutta-Nyström method which can be applied directly to second order ODEs.
In particular, this method allows the acceleration equation to depend on the velocity.

## References

```
@article{fine1987low,
  title={Low order practical {R}unge-{K}utta-{N}ystr{\"o}m methods},
  author={Fine, Jerry Michael},
  journal={Computing},
  volume={38},
  number={4},
  pages={281--297},
  year={1987},
  publisher={Springer}
}
```
"""
struct FineRKN5 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    Nystrom4VelocityIdependent

A 4th order explicit Runkge-Kutta-Nyström method. Used directly on second order ODEs, where the acceleration is independent from velocity (ODE Problem is not dependent on the first derivative).

More efficient then [`Nystrom4`](@ref) on velocity independent problems, since less evaluations are needed.

Fixed time steps only.

## References

E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
Springer-Verlag.
"""
struct Nystrom4VelocityIndependent <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    IRKN4

Improves Runge-Kutta-Nyström method of order four, which minimizes the amount of evaluated functions in each step. Fixed time steps only.

Second order ODE should not be dependent on the first derivative.

Recommended for smooth problems with expensive functions to evaluate.

## References

@article{rabiei2012numerical,
title={Numerical Solution of Second-Order Ordinary Differential Equations by Improved Runge-Kutta Nystrom Method},
author={Rabiei, Faranak and Ismail, Fudziah and Norazak, S and Emadi, Saeid},
publisher={Citeseer}
}
"""
struct IRKN4 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    Nystrom5VelocityIndependent

A 5th order explicit Runkge-Kutta-Nyström method. Used directly on second order ODEs, where the acceleration is independent from velocity (ODE Problem is not dependent on the first derivative).
Fixed time steps only.

## References

E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
Springer-Verlag.
"""
struct Nystrom5VelocityIndependent <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    DPRKN4

4th order explicit Runge-Kutta-Nyström methods. The second order ODE should not depend on the first derivative.

## References

@article{Dormand1987FamiliesOR,
title={Families of Runge-Kutta-Nystrom Formulae},
author={J. R. Dormand and Moawwad E. A. El-Mikkawy and P. J. Prince},
journal={Ima Journal of Numerical Analysis},
year={1987},
volume={7},
pages={235-250}
}
"""
struct DPRKN4 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN5

5th order explicit Runge-Kutta-Nyström mehod. The second order ODE should not depend on the first derivative.

## References

@article{Bettis1973ARN,
title={A Runge-Kutta Nystrom algorithm},
author={Dale G. Bettis},
journal={Celestial mechanics},
year={1973},
volume={8},
pages={229-233},
publisher={Springer}
}
"""
struct DPRKN5 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN6

6th order explicit Runge-Kutta-Nyström method. The second order ODE should not depend on the first derivative. Free 6th order interpolant.

## References

@article{dormand1987runge,
title={Runge-kutta-nystrom triples},
author={Dormand, JR and Prince, PJ},
journal={Computers \\& Mathematics with Applications},
volume={13},
number={12},
pages={937--949},
year={1987},
publisher={Elsevier}
}
"""
struct DPRKN6 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN6FM

6th order explicit Runge-Kutta-Nyström method. The second order ODE should not depend on the first derivative.

Compared to [`DPRKN6`](@ref), this method has smaller truncation error coefficients which leads to performance gain
when only the main solution points are considered.

## References

@article{Dormand1987FamiliesOR,
title={Families of Runge-Kutta-Nystrom Formulae},
author={J. R. Dormand and Moawwad E. A. El-Mikkawy and P. J. Prince},
journal={Ima Journal of Numerical Analysis},
year={1987},
volume={7},
pages={235-250}
}
"""
struct DPRKN6FM <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN8

8th order explicit Runge-Kutta-Nyström method. The second order ODE should not depend on the first derivative.

Not as efficient as [`DPRKN12`](@ref) when high accuracy is needed, however this solver is competitive with
[`DPRKN6`](@ref) at lax tolerances and, depending on the problem, might be a good option between performance and accuracy.

## References

@article{dormand1987high,
title={High-order embedded Runge-Kutta-Nystrom formulae},
author={Dormand, JR and El-Mikkawy, MEA and Prince, PJ},
journal={IMA Journal of Numerical Analysis},
volume={7},
number={4},
pages={423--430},
year={1987},
publisher={Oxford University Press}
}
"""
struct DPRKN8 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN12

12th order explicit Rugne-Kutta-Nyström method. The second order ODE should not depend on the first derivative.

Most efficient when high accuracy is needed.

## References

@article{dormand1987high,
title={High-order embedded Runge-Kutta-Nystrom formulae},
author={Dormand, JR and El-Mikkawy, MEA and Prince, PJ},
journal={IMA Journal of Numerical Analysis},
volume={7},
number={4},
pages={423--430},
year={1987},
publisher={Oxford University Press}
}
"""
struct DPRKN12 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    ERKN4

Embedded 4(3) pair of explicit Runge-Kutta-Nyström methods. Integrates the periodic properties of the harmonic oscillator exactly.

The second order ODE should not depend on the first derivative.

Uses adaptive step size control. This method is extra efficient on periodic problems.

## References

@article{demba2017embedded,
title={An Embedded 4 (3) Pair of Explicit Trigonometrically-Fitted Runge-Kutta-Nystr{\"o}m Method for Solving Periodic Initial Value Problems},
author={Demba, MA and Senu, N and Ismail, F},
journal={Applied Mathematical Sciences},
volume={11},
number={17},
pages={819--838},
year={2017}
}
"""
struct ERKN4 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    ERKN5

Embedded 5(4) pair of explicit Runge-Kutta-Nyström methods. Integrates the periodic properties of the harmonic oscillator exactly.

The second order ODE should not depend on the first derivative.

Uses adaptive step size control. This method is extra efficient on periodic problems.

## References

@article{demba20165,
title={A 5 (4) Embedded Pair of Explicit Trigonometrically-Fitted Runge--Kutta--Nystr{\"o}m Methods for the Numerical Solution of Oscillatory Initial Value Problems},
author={Demba, Musa A and Senu, Norazak and Ismail, Fudziah},
journal={Mathematical and Computational Applications},
volume={21},
number={4},
pages={46},
year={2016},
publisher={Multidisciplinary Digital Publishing Institute}
}
"""
struct ERKN5 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    ERKN7

Embedded pair of explicit Runge-Kutta-Nyström methods. Integrates the periodic properties of the harmonic oscillator exactly.

The second order ODE should not depend on the first derivative.

Uses adaptive step size control. This method is extra efficient on periodic Problems.

## References

@article{SimosOnHO,
title={On high order Runge-Kutta-Nystr{\"o}m pairs},
author={Theodore E. Simos and Ch. Tsitouras},
journal={J. Comput. Appl. Math.},
volume={400},
pages={113753}
}
"""
struct ERKN7 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
3 stage fourth order Runge-Kutta Nystrom method to solve second order linear inhomogeneous IVPs.

Does not include an adaptive method. Solves for for d-dimensional differential systems of second order linear inhomogeneous equations.

!!! warn

    This method is only fourth order for these systems, the method is second order otherwise!

## References

@article{MONTIJANO2024115533,
title = {Explicit Runge–Kutta–Nyström methods for the numerical solution of second order linear inhomogeneous IVPs},
author = {J.I. Montijano and L. Rández and M. Calvo},
journal = {Journal of Computational and Applied Mathematics},
volume = {438},
pages = {115533},
year = {2024},
}
"""
struct RKN4 <: OrdinaryDiffEqAlgorithm end