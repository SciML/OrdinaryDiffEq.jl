alg_order(alg::KYK2014DGSSPRK_3S2) = 2
alg_order(alg::ORK256) = 2
alg_order(alg::CarpenterKennedy2N54) = 4
alg_order(alg::SSPRK53_2N1) = 3
alg_order(alg::SSPRK53_2N2) = 3
alg_order(alg::NDBLSRK124) = 4
alg_order(alg::NDBLSRK134) = 4
alg_order(alg::NDBLSRK144) = 4
alg_order(alg::CFRLDDRK64) = 4
alg_order(alg::DGLDDRK73_C) = 3
alg_order(alg::TSLDDRK74) = 4
alg_order(alg::DGLDDRK84_C) = 4
alg_order(alg::DGLDDRK84_F) = 4
alg_order(alg::SHLDDRK64) = 4
alg_order(alg::RK46NL) = 4
alg_order(alg::ParsaniKetchesonDeconinck3S32) = 2
alg_order(alg::ParsaniKetchesonDeconinck3S82) = 2
alg_order(alg::ParsaniKetchesonDeconinck3S53) = 3
alg_order(alg::ParsaniKetchesonDeconinck3S173) = 3
alg_order(alg::ParsaniKetchesonDeconinck3S94) = 4
alg_order(alg::ParsaniKetchesonDeconinck3S184) = 4
alg_order(alg::ParsaniKetchesonDeconinck3S105) = 5
alg_order(alg::ParsaniKetchesonDeconinck3S205) = 5
alg_order(alg::CKLLSRK43_2) = 3
alg_order(alg::CKLLSRK54_3C) = 4
alg_order(alg::CKLLSRK95_4S) = 5
alg_order(alg::CKLLSRK95_4C) = 5
alg_order(alg::CKLLSRK95_4M) = 5
alg_order(alg::CKLLSRK54_3C_3R) = 4
alg_order(alg::CKLLSRK54_3M_3R) = 4
alg_order(alg::CKLLSRK54_3N_3R) = 4
alg_order(alg::CKLLSRK85_4C_3R) = 5
alg_order(alg::CKLLSRK85_4M_3R) = 5
alg_order(alg::RDPK3Sp35) = 3
alg_order(alg::RDPK3SpFSAL35) = 3
alg_order(alg::RDPK3Sp49) = 4
alg_order(alg::RDPK3SpFSAL49) = 4
alg_order(alg::RDPK3Sp510) = 5
alg_order(alg::RDPK3SpFSAL510) = 5
alg_order(alg::CKLLSRK85_4P_3R) = 5
alg_order(alg::CKLLSRK85_4FM_4R) = 5
alg_order(alg::CKLLSRK54_3N_4R) = 4
alg_order(alg::CKLLSRK75_4M_5R) = 5
alg_order(alg::CKLLSRK54_3M_4R) = 4
alg_order(alg::CKLLSRK65_4M_4R) = 5

alg_order(alg::SSPRK22) = 2
alg_order(alg::SSPRKMSVS32) = 2
alg_order(alg::SSPRK33) = 3
alg_order(alg::KYKSSPRK42) = 2
alg_order(alg::SSPRK53) = 3
alg_order(alg::SSPRK53_H) = 3
alg_order(alg::SSPRK63) = 3
alg_order(alg::SSPRK73) = 3
alg_order(alg::SSPRK83) = 3
alg_order(alg::SSPRK43) = 3
alg_order(alg::SSPRK432) = 3
alg_order(alg::SSPRKMSVS43) = 3
alg_order(alg::SSPRK932) = 3
alg_order(alg::SSPRK54) = 4
alg_order(alg::SSPRK104) = 4
alg_order(alg::SHLDDRK52) = 2
alg_order(alg::SHLDDRK_2N) = 4

isfsal(alg::ORK256) = false
isfsal(alg::CarpenterKennedy2N54) = false
isfsal(alg::SSPRK53_2N1) = false
isfsal(alg::SSPRK53_2N2) = false
isfsal(alg::DGLDDRK84_F) = false
isfsal(alg::DGLDDRK73_C) = false
isfsal(alg::DGLDDRK84_C) = false
isfsal(alg::RDPK3Sp35) = false
isfsal(alg::RDPK3Sp49) = false
isfsal(alg::RDPK3Sp510) = false
isfsal(alg::SHLDDRK64) = false
isfsal(alg::NDBLSRK134) = false
isfsal(alg::NDBLSRK124) = false
isfsal(alg::NDBLSRK144) = false

isfsal(alg::SSPRK22) = false
isfsal(alg::SSPRK33) = false
isfsal(alg::SSPRK53) = false
isfsal(alg::SSPRK53_H) = false
isfsal(alg::SSPRK63) = false
isfsal(alg::SSPRK73) = false
isfsal(alg::SSPRK83) = false
isfsal(alg::SSPRK43) = false
isfsal(alg::SSPRK432) = false
isfsal(alg::SSPRK932) = false
isfsal(alg::SSPRK54) = false
isfsal(alg::SSPRK104) = false

uses_uprev(alg::ORK256, adaptive::Bool) = false
uses_uprev(alg::SHLDDRK64, adaptive::Bool) = false
uses_uprev(alg::CarpenterKennedy2N54, adaptive::Bool) = false
uses_uprev(alg::NDBLSRK124, adaptive::Bool) = false
uses_uprev(alg::NDBLSRK134, adaptive::Bool) = false
uses_uprev(alg::DGLDDRK84_F, adaptive::Bool) = false
uses_uprev(alg::NDBLSRK144, adaptive::Bool) = false
uses_uprev(alg::DGLDDRK84_C, adaptive::Bool) = false
uses_uprev(alg::TSLDDRK74, adaptive::Bool) = false
uses_uprev(alg::CFRLDDRK64, adaptive::Bool) = false
uses_uprev(alg::DGLDDRK73_C, adaptive::Bool) = false
uses_uprev(alg::CKLLSRK43_2, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3C, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK95_4S, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK95_4C, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK95_4M, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3C_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3M_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3N_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK85_4C_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK85_4M_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK85_4P_3R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3N_4R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK54_3M_4R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK65_4M_4R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK85_4FM_4R, adaptive::Bool) = adaptive
uses_uprev(alg::CKLLSRK75_4M_5R, adaptive::Bool) = adaptive

"""
    ssp_coefficient(alg)

Return the SSP coefficient of the ODE algorithm `alg`. If one time step of size
`dt` with `alg` can be written as a convex combination of explicit Euler steps
with step sizes `cᵢ * dt`, the SSP coefficient is the minimal value of `1/cᵢ`.

# Examples

```julia-repl
julia> ssp_coefficient(SSPRK104())
6
```
"""

ssp_coefficient(alg::SSPRK53_2N1) = 2.18
ssp_coefficient(alg::SSPRK53_2N2) = 2.148
ssp_coefficient(alg::SSPRK53) = 2.65
ssp_coefficient(alg::SSPRK53_H) = 2.65
ssp_coefficient(alg::SSPRK63) = 3.518
ssp_coefficient(alg::SSPRK73) = 4.2879
ssp_coefficient(alg::SSPRK83) = 5.107
ssp_coefficient(alg::SSPRK43) = 2
ssp_coefficient(alg::SSPRK432) = 2
ssp_coefficient(alg::KYKSSPRK42) = 2.459
ssp_coefficient(alg::SSPRK932) = 6
ssp_coefficient(alg::SSPRK54) = 1.508
ssp_coefficient(alg::SSPRK104) = 6
ssp_coefficient(alg::KYK2014DGSSPRK_3S2) = 0.8417
ssp_coefficient(alg::SSPRK33) = 1
ssp_coefficient(alg::SSPRK22) = 1