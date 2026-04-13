using StochasticDiffEq, Test

SRIW1_tab = constructSRIW1()
@test minimum(checkSRIOrder(SRIW1_tab))

SRA1_tab = constructSRA1()
@test minimum(checkSRAOrder(SRA1_tab))

DRI1_tab = constructDRI1()
@test minimum(checkRIOrder(DRI1_tab))

RI1_tab = constructRI1()
@test minimum(checkRIOrder(RI1_tab))

RI3_tab = constructRI3()
@test minimum(checkRIOrder(RI3_tab))

RI5_tab = constructRI5()
@test minimum(checkRIOrder(RI5_tab))

RI6_tab = constructRI6()
@test minimum(checkRIOrder(RI6_tab))

RDI1WM_tab = constructRDI1WM()
@test minimum(checkRIOrder(RDI1WM_tab, ps = 1))

RDI2WM_tab = constructRDI2WM()
@test minimum(checkRIOrder(RDI2WM_tab))

RDI3WM_tab = constructRDI3WM()
@test minimum(checkRIOrder(RDI3WM_tab))

RDI4WM_tab = constructRDI4WM()
@test minimum(checkRIOrder(RDI4WM_tab))

RS1_tab = constructRS1()
@test minimum(checkRSOrder(RS1_tab))

RS2_tab = constructRS2()
@test minimum(checkRSOrder(RS2_tab))

# Komori's NON
NON_tab = constructNON()
@test minimum(checkNONOrder(NON_tab))

NON2_tab = constructNON2()
@test minimum(checkNONOrder(NON2_tab))
