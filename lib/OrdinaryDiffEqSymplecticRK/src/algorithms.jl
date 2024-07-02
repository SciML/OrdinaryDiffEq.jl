struct SymplecticEuler <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{verlet1967computer,
title={Computer" experiments" on classical fluids. I. Thermodynamical properties of Lennard-Jones molecules},
author={Verlet, Loup},
journal={Physical review},
volume={159},
number={1},
pages={98},
year={1967},
publisher={APS}
}
"""

struct VelocityVerlet <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{verlet1967computer,
title={Computer" experiments" on classical fluids. I. Thermodynamical properties of Lennard-Jones molecules},
author={Verlet, Loup},
journal={Physical review},
volume={159},
number={1},
pages={98},
year={1967},
publisher={APS}
}
"""

struct VerletLeapfrog <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{verlet1967computer,
title={Computer" experiments" on classical fluids. I. Thermodynamical properties of Lennard-Jones molecules},
author={Verlet, Loup},
journal={Physical review},
volume={159},
number={1},
pages={98},
year={1967},
publisher={APS}
}
"""

struct PseudoVerletLeapfrog <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1992accuracy,
title={The accuracy of symplectic integrators},
author={McLachlan, Robert I and Atela, Pau},
journal={Nonlinearity},
volume={5},
number={2},
pages={541},
year={1992},
publisher={IOP Publishing}
}
"""

struct McAte2 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{ruth1983canonical,
title={A canonical integration technique},
author={Ruth, Ronald D},
journal={IEEE Trans. Nucl. Sci.},
volume={30},
number={CERN-LEP-TH-83-14},
pages={2669--2671},
year={1983}
}
"""

struct Ruth3 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1992accuracy,
title={The accuracy of symplectic integrators},
author={McLachlan, Robert I and Atela, Pau},
journal={Nonlinearity},
volume={5},
number={2},
pages={541},
year={1992},
publisher={IOP Publishing}
}
"""

struct McAte3 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{candy1991symplectic,
title={A symplectic integration algorithm for separable Hamiltonian functions},
author={Candy, J and Rozmus, W},
journal={Journal of Computational Physics},
volume={92},
number={1},
pages={230--256},
year={1991},
publisher={Elsevier}
}
"""

struct CandyRoz4 <: OrdinaryDiffEqPartitionedAlgorithm end

struct McAte4 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{sanz1993symplectic,
title={Symplectic numerical methods for Hamiltonian problems},
author={Sanz-Serna, Jes{\'u}s Maria and Calvo, Mari-Paz},
journal={International Journal of Modern Physics C},
volume={4},
number={02},
pages={385--392},
year={1993},
publisher={World Scientific}
}
"""

struct CalvoSanz4 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1992accuracy,
title={The accuracy of symplectic integrators},
author={McLachlan, Robert I and Atela, Pau},
journal={Nonlinearity},
volume={5},
number={2},
pages={541},
year={1992},
publisher={IOP Publishing}
}
"""

struct McAte42 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1992accuracy,
title={The accuracy of symplectic integrators},
author={McLachlan, Robert I and Atela, Pau},
journal={Nonlinearity},
volume={5},
number={2},
pages={541},
year={1992},
publisher={IOP Publishing}
}
"""

struct McAte5 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{yoshida1990construction,
title={Construction of higher order symplectic integrators},
author={Yoshida, Haruo},
journal={Physics letters A},
volume={150},
number={5-7},
pages={262--268},
year={1990},
publisher={Elsevier}
}
"""

struct Yoshida6 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{kahan1997composition,
title={Composition constants for raising the orders of unconventional schemes for ordinary differential equations},
author={Kahan, William and Li, Ren-Cang},
journal={Mathematics of computation},
volume={66},
number={219},
pages={1089--1099},
year={1997}
}
"""

struct KahanLi6 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1995numerical,
title={On the numerical integration of ordinary differential equations by symmetric composition methods},
author={McLachlan, Robert I},
journal={SIAM Journal on Scientific Computing},
volume={16},
number={1},
pages={151--168},
year={1995},
publisher={SIAM}
}
"""

struct McAte8 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{kahan1997composition,
title={Composition constants for raising the orders of unconventional schemes for ordinary differential equations},
author={Kahan, William and Li, Ren-Cang},
journal={Mathematics of computation},
volume={66},
number={219},
pages={1089--1099},
year={1997}
}
"""

struct KahanLi8 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{sofroniou2005derivation,
title={Derivation of symmetric composition constants for symmetric integrators},
author={Sofroniou, Mark and Spaletta, Giulia},
journal={Optimization Methods and Software},
volume={20},
number={4-5},
pages={597--613},
year={2005},
publisher={Taylor \\& Francis}
}
"""
struct SofSpa10 <: OrdinaryDiffEqPartitionedAlgorithm end