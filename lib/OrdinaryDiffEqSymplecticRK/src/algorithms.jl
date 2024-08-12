@doc generic_solver_docstring("DESCRIPTION TBD",
        "SymplecticEuler",
        "Symplectic Runge-Kutta Methods",
        "REF TBD", "", "")
struct SymplecticEuler <: OrdinaryDiffEqPartitionedAlgorithm end

verlet1967 = """
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

@doc generic_solver_docstring("DESCRIPTION TBD",
        "VelocityVerlet",
        "Symplectic Runge-Kutta Methods",
        verlet1967, "", "")
struct VelocityVerlet <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
        "VerletLeapfrog",
        "Symplectic Runge-Kutta Methods",
        verlet1967, "", "")
struct VerletLeapfrog <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
        "PseudoVerletLeapfrog",
        "Symplectic Runge-Kutta Methods",
        verlet1967, "", "")
struct PseudoVerletLeapfrog <: OrdinaryDiffEqPartitionedAlgorithm end

mclachlan1992 = """
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

@doc generic_solver_docstring("DESCRIPTION TBD",
        "McAte2",
        "Symplectic Runge-Kutta Methods",
        mclachlan1992, "", "")
struct McAte2 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
    "Ruth3",
    "Symplectic Runge-Kutta Methods",
    """@article{ruth1983canonical,
    title={A canonical integration technique},
    author={Ruth, Ronald D},
    journal={IEEE Trans. Nucl. Sci.},
    volume={30},
    number={CERN-LEP-TH-83-14},
    pages={2669--2671},
    year={1983}}""", "", "")
struct Ruth3 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
        "McAte3",
        "Symplectic Runge-Kutta Methods",
        mclachlan1992, "", "")
struct McAte3 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
    "CandyRoz4",
    "Symplectic Runge-Kutta Methods",
    """@article{candy1991symplectic,
    itle={A symplectic integration algorithm for separable Hamiltonian functions},
    uthor={Candy, J and Rozmus, W},
    ournal={Journal of Computational Physics},
    olume={92},
    umber={1},
    ages={230--256},
    ear={1991},
    ublisher={Elsevier}}""", "", "")
struct CandyRoz4 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
        "McAte4",
        "Symplectic Runge-Kutta Methods",
        mclachlan1992, "", "")
struct McAte4 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
    "CalvoSanz4",
    "Symplectic Runge-Kutta Methods",
    """@article{sanz1993symplectic,
    title={Symplectic numerical methods for Hamiltonian problems},
    author={Sanz-Serna, Jes{\'u}s Maria and Calvo, Mari-Paz},
    journal={International Journal of Modern Physics C},
    volume={4},
    number={02},
    pages={385--392},
    year={1993},
    publisher={World Scientific}
    }""", "", "")
struct CalvoSanz4 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
        "McAte42",
        "Symplectic Runge-Kutta Methods",
        mclachlan1992, "", "")
struct McAte42 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
        "McAte5",
        "Symplectic Runge-Kutta Methods",
        mclachlan1992, "", "")
struct McAte5 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
    "Yoshida6",
    "Symplectic Runge-Kutta Methods",
    """@article{yoshida1990construction,
    title={Construction of higher order symplectic integrators},
    author={Yoshida, Haruo},
    journal={Physics letters A},
    volume={150},
    number={5-7},
    pages={262--268},
    year={1990},
    publisher={Elsevier}}""", "", "")
struct Yoshida6 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
    "KahanLi6",
    "Symplectic Runge-Kutta Methods",
    """@article{yoshida1990construction,
    title={Construction of higher order symplectic integrators},
    author={Yoshida, Haruo},
    journal={Physics letters A},
    volume={150},
    number={5-7},
    pages={262--268},
    year={1990},
    publisher={Elsevier}}""", "", "")
struct KahanLi6 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
    "McAte8",
    "Symplectic Runge-Kutta Methods",
    """@article{mclachlan1995numerical,
    title={On the numerical integration of ordinary differential equations by symmetric composition methods},
    author={McLachlan, Robert I},
    journal={SIAM Journal on Scientific Computing},
    volume={16},
    number={1},
    pages={151--168},
    year={1995},
    publisher={SIAM}
    }""", "", "")
struct McAte8 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
    "KahanLi8",
    "Symplectic Runge-Kutta Methods",
    """@article{kahan1997composition,
    title={Composition constants for raising the orders of unconventional schemes for ordinary differential equations},
    author={Kahan, William and Li, Ren-Cang},
    journal={Mathematics of computation},
    volume={66},
    number={219},
    pages={1089--1099},
    year={1997}}""", "", "")
struct KahanLi8 <: OrdinaryDiffEqPartitionedAlgorithm end

@doc generic_solver_docstring("DESCRIPTION TBD",
    "SofSpa10",
    "Symplectic Runge-Kutta Methods",
    """@article{sofroniou2005derivation,
    title={Derivation of symmetric composition constants for symmetric integrators},
    author={Sofroniou, Mark and Spaletta, Giulia},
    journal={Optimization Methods and Software},
    volume={20},
    number={4-5},
    pages={597--613},
    year={2005},
    publisher={Taylor \\& Francis}}""", "", "")
struct SofSpa10 <: OrdinaryDiffEqPartitionedAlgorithm end
