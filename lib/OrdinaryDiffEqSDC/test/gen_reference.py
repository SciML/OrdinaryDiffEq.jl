#!/usr/bin/env python3
"""Generate golden reference coefficients + expected SDC orders from qmat."""
import json
import sys

sys.path.insert(0, "qmat")

import numpy as np
from qmat.qcoeff.collocation import Collocation
from qmat import genQDeltaCoeffs
from qmat.utils.sdc import getOrderSDC, errorDahlquistSDC
from qmat.utils.num import numericalOrder

NODE_TYPES = ["LEGENDRE", "EQUID"]
QUAD_TYPES = ["GAUSS", "RADAU-LEFT", "RADAU-RIGHT", "LOBATTO"]
SWEEPERS = ["BE", "FE", "TRAP", "LU", "PIC", "BEPAR", "MIN-SR-NS"]

out = {"coeffs": [], "orders": [], "measured": []}

for nodeType in NODE_TYPES:
    for quadType in QUAD_TYPES:
        for M in [2, 3, 4, 5]:
            coll = Collocation(M, nodeType, quadType)
            nodes, weights, Q = coll.genCoeffs()
            entry = {
                "M": M,
                "nodeType": nodeType,
                "quadType": quadType,
                "collOrder": coll.order,
                "nodes": nodes.tolist(),
                "weights": weights.tolist(),
                "Q": Q.tolist(),
                "QDelta": {},
            }
            for s in SWEEPERS:
                try:
                    QD = genQDeltaCoeffs(s, qGen=coll, nodes=nodes, Q=Q)
                except Exception as e:  # noqa: BLE001
                    entry["QDelta"][s] = f"ERROR: {e}"
                    continue
                entry["QDelta"][s] = np.asarray(QD).tolist()
            out["coeffs"].append(entry)

# Expected orders (qmat's own predictor) for the combinations we will test.
lam, u0, T = 1j, 1, 2 * np.pi


def nStepsForTest(order):
    nSteps = [1, 2, 4]
    if order == 1:
        nSteps = [64, 128, 256]
    elif order == 2:
        nSteps = [32, 64, 128]
    elif order == 3:
        nSteps = [16, 32, 64]
    elif order in [4, 5]:
        nSteps = [8, 16, 32]
    elif order in [6, 7]:
        nSteps = [4, 8, 16]
    return nSteps


for nodeType in NODE_TYPES:
    for quadType in QUAD_TYPES:
        for M in [2, 3, 4]:
            coll = Collocation(M, nodeType, quadType)
            nodes, weights, Q = coll.genCoeffs()
            for s in SWEEPERS:
                QD = np.asarray(genQDeltaCoeffs(s, qGen=coll, nodes=nodes, Q=Q))
                for K in [1, 2, 3, 4]:
                    for su in ["QUADRATURE", "LASTNODE"]:
                        if su == "LASTNODE" and quadType not in (
                            "RADAU-RIGHT",
                            "LOBATTO",
                        ):
                            continue
                        p = getOrderSDC(coll, K, s, su)
                        nSteps = nStepsForTest(p)
                        w = weights if su == "QUADRATURE" else None
                        err = [
                            errorDahlquistSDC(lam, u0, T, nS, K, Q, QD, w)
                            for nS in nSteps
                        ]
                        obs, rmse = numericalOrder(nSteps, err)
                        out["orders"].append(
                            {
                                "M": M,
                                "nodeType": nodeType,
                                "quadType": quadType,
                                "sweeper": s,
                                "K": K,
                                "stepUpdate": su,
                                "predicted": p,
                                "nSteps": nSteps,
                                "err": [float(e) for e in err],
                                "observed": float(obs),
                                "rmse": float(rmse),
                            }
                        )

with open("qmat_reference.json", "w") as fh:
    json.dump(out, fh, indent=1)
print("coeff entries:", len(out["coeffs"]))
print("order entries:", len(out["orders"]))
