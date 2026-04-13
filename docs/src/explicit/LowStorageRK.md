```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqLowStorageRK

Low-storage Runge-Kutta methods are specialized explicit schemes designed to minimize memory requirements while maintaining high-order accuracy. These methods are essential for large-scale computational fluid dynamics and wave propagation problems where memory constraints are critical.

## Key Properties

Low-storage RK methods provide:

  - **Drastically reduced memory usage** (typically 2-4 registers vs 7-10 for standard RK)
  - **High-order accuracy** comparable to standard RK methods
  - **Preservation of important stability properties** (low dissipation/dispersion)
  - **Scalability** to very large PDE discretizations

## When to Use Low-Storage RK Methods

These methods are recommended for:

  - **Large-scale PDE discretizations** where memory is the limiting factor
  - **Computational fluid dynamics** and wave propagation simulations
  - **High-performance computing** applications with memory constraints
  - **GPU computations** where memory bandwidth is critical
  - **Compressible flow simulations** and aerodynamics
  - **Seismic wave propagation** and acoustic simulations
  - **Problems with millions or billions of unknowns**

## Memory Efficiency Comparison

**Registers** refer to the number of copies of the `u0` vector that must be stored in memory during integration:

  - **Standard Tsit5**: ~9 registers (copies of the state vector)
  - **Low-storage methods**: 2-4 registers (copies of the state vector)
  - **Practical example**: If `u0` is from a PDE semi-discretization requiring 2 GB, then Tsit5 needs 18 GB of working memory, while a 2-register method only needs 4 GB and can achieve the same order
  - **Trade-off**: These methods achieve memory reduction by being less computationally efficient, trading compute performance for lower memory requirements

## Solver Selection Guide

### General-purpose low-storage

  - **`CarpenterKennedy2N54`**: Fourth-order, 5-stage, excellent general choice
  - **`RDPK3Sp510`**: Fifth-order with only 3 registers, very memory efficient

### Wave propagation optimized

  - **`ORK256`**: Second-order, 5-stage, optimized for wave equations
  - **`CFRLDDRK64`**: Low-dissipation and low-dispersion variant
  - **`TSLDDRK74`**: Seventh-order for high accuracy wave propagation

### Discontinuous Galerkin optimized

  - **`DGLDDRK73_C`**: Optimized for DG discretizations (constrained)
  - **`DGLDDRK84_C`**, **`DGLDDRK84_F`**: Fourth-order DG variants

### Specialized high-order

  - **`NDBLSRK124`**, **`NDBLSRK144`**: Multi-stage fourth-order methods
  - **`SHLDDRK64`**: Low dissipation and dispersion properties
  - **`RK46NL`**: Six-stage fourth-order method

### Computational fluid dynamics

  - **Carpenter-Kennedy-Lewis series** (`CKLLSRK*`): Optimized for Navier-Stokes equations
  - **Parsani-Ketcheson-Deconinck series** (`ParsaniKetcheson*`): CFD-optimized variants
  - **Ranocha-Dalcin-Parsani-Ketcheson series** (`RDPK*`): Modern CFD methods

## Performance Considerations

  - **Use only when memory-bound**: Standard RK methods are often more efficient when memory is not limiting
  - **Best for large systems**: Most beneficial for problems with >10⁶ unknowns
  - **GPU acceleration**: Particularly effective on memory-bandwidth limited hardware

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqLowStorageRK", "CarpenterKennedy2N54")
```

## Full list of solvers

```@docs
ORK256
DGLDDRK73_C
CarpenterKennedy2N54
NDBLSRK124
NDBLSRK144
CFRLDDRK64
TSLDDRK74
DGLDDRK84_C
DGLDDRK84_F
SHLDDRK64
RK46NL
ParsaniKetchesonDeconinck3S32
ParsaniKetchesonDeconinck3S82
ParsaniKetchesonDeconinck3S53
ParsaniKetchesonDeconinck3S173
ParsaniKetchesonDeconinck3S94
ParsaniKetchesonDeconinck3S184
ParsaniKetchesonDeconinck3S105
ParsaniKetchesonDeconinck3S205
CKLLSRK43_2
CKLLSRK54_3C
CKLLSRK95_4S
CKLLSRK95_4C
CKLLSRK95_4M
CKLLSRK54_3C_3R
CKLLSRK54_3M_3R
CKLLSRK54_3N_3R
CKLLSRK85_4C_3R
CKLLSRK85_4M_3R
CKLLSRK85_4P_3R
CKLLSRK54_3N_4R
CKLLSRK54_3M_4R
CKLLSRK65_4M_4R
CKLLSRK85_4FM_4R
CKLLSRK75_4M_5R
RDPK3Sp35
RDPK3SpFSAL35
RDPK3Sp49
RDPK3SpFSAL49
RDPK3Sp510
RDPK3SpFSAL510
HSLDDRK64
NDBLSRK134
SHLDDRK_2N
SHLDDRK52
```
