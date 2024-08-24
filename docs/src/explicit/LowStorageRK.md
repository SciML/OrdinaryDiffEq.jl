```@meta
CollapsedDocStrings = true
```
# OrdinaryDiffEqLowStorageRK

These methods are designed to have reduced register requirements, allowing for larger sets of ODEs to more
easily fit into RAM. For example, while the 5th order Tsit5 requires around 9 concurrent instantiations of the
ODE state `u`, `RDPK3Sp510` can achieve 5th order with 3 registers, effectively requiring 1/3 of the memory.
However, there are some efficiency trade-offs used in the design of the low-storage RK methods, and thus they
are generally only recommended in situations which are RAM bound, like large-scale PDE discretizations.

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
KYK2014DGSSPRK_3S2
```
