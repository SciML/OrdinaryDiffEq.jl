# Stability Limit Detection (STALD) for BDF methods
#
# This implements the algorithm from CVODE (Hindmarsh 1992, 1995) for detecting
# when the BDF method at orders 3-5 is operating near its stability boundary.
# BDF orders 1-2 are A-stable (unconditionally stable for all eigenvalues in the
# left half-plane). Orders 3-5 are only alpha-stable, meaning their stability
# regions exclude portions of the left half-plane near the imaginary axis.
#
# When the dominant eigenvalue's step-size-scaled magnitude approaches the stability
# boundary, the scaled derivative norms exhibit a characteristic geometric pattern.
# This algorithm detects that pattern and forces order reduction when needed.

"""
    StabilityLimitDetectionState{T}

State for the CVODE-style Stability Limit Detection (STALD) algorithm.
Tracks scaled derivative data over a sliding window of 5 steps to detect
when BDF orders 3-5 are near their stability boundaries.

Parameterized on element type `T` (typically `real(eltype(u))`) to avoid
hardcoded Float64 and support other floating-point types.
"""
mutable struct StabilityLimitDetectionState{T}
    ssdat::Matrix{T}    # 5×3 matrix: ssdat[i,k] stores squared scaled derivative norms
    # i = step index (1=newest, 5=oldest), k = derivative level (1=q-2, 2=q-1, 3=q)
    nscon::Int           # consecutive steps at same order
    nor::Int             # total number of STALD-triggered order reductions
    enabled::Bool        # whether STALD is active
    last_order::Int      # order used on the previous step
    # Tuning constants (matching CVODE defaults)
    rrcut::T             # cutoff for characteristic root magnitude
    vrrtol::T            # tolerance for variance of ratios
    vrrt2::T             # secondary variance tolerance
    sqtol::T             # tolerance for quartic residual
    rrtol::T             # tolerance for rr cross-verification
    tiny::T              # tiny value to avoid division by zero
    # Pre-allocated workspace for _stald_detect (avoids MVector/MMatrix allocations)
    smax::Vector{T}      # length 3
    ssmax::Vector{T}     # length 3
    rav::Vector{T}       # length 3
    vrat::Vector{T}      # length 3
    rat::Matrix{T}       # 4×3
    qc::Matrix{T}        # 5×3
    qco::Matrix{T}       # 5×3
    drr::Vector{T}       # length 3
    rrc::Vector{T}       # length 3
    sigsq::Vector{T}     # length 3
end

function StabilityLimitDetectionState(
        ::Type{T} = Float64;
        enabled::Bool = true,
        rrcut = T(0.98),
        vrrtol = T(1.0e-4),
        vrrt2 = T(5.0e-4),
        sqtol = T(1.0e-3),
        rrtol = T(1.0e-2),
        tiny = T(1.0e-90),
    ) where {T}
    return StabilityLimitDetectionState{T}(
        zeros(T, 5, 3),
        0,
        0,
        enabled,
        0,
        T(rrcut),
        T(vrrtol),
        T(vrrt2),
        T(sqtol),
        T(rrtol),
        T(tiny),
        zeros(T, 3),       # smax
        zeros(T, 3),       # ssmax
        zeros(T, 3),       # rav
        zeros(T, 3),       # vrat
        zeros(T, 4, 3),    # rat
        zeros(T, 5, 3),    # qc
        zeros(T, 5, 3),    # qco
        zeros(T, 3),       # drr
        zeros(T, 3),       # rrc
        zeros(T, 3),       # sigsq
    )
end

"""
    stald_collect_data!(stald, order, terkm2, terkm1, terk)

Collect scaled derivative data for STALD after a successful step.
Only collects data when order >= 3 (orders 1-2 are A-stable).

Arguments:
- `stald`: StabilityLimitDetectionState
- `order`: current BDF order
- `terkm2, terkm1, terk`: weighted norm estimates of h^(k-2)*y^(k-2), h^(k-1)*y^(k-1), h^k*y^k
"""
function stald_collect_data!(
        stald::StabilityLimitDetectionState, order::Int,
        terkm2::Real, terkm1::Real, terk::Real
    )
    if !stald.enabled || order < 3
        return nothing
    end

    # Shift old data down (newest at index 1)
    for i in 5:-1:2
        for k in 1:3
            stald.ssdat[i, k] = stald.ssdat[i - 1, k]
        end
    end

    # Store squared norms at the three derivative levels
    # Using the same convention as CVODE:
    # k=1 corresponds to order q-2 estimate
    # k=2 corresponds to order q-1 estimate
    # k=3 corresponds to order q estimate
    stald.ssdat[1, 1] = terkm2 * terkm2
    stald.ssdat[1, 2] = terkm1 * terkm1
    stald.ssdat[1, 3] = terk * terk

    return nothing
end

"""
    stald_check!(stald, order) -> Bool

Check for stability limit violation and return true if order should be reduced.
Called after each successful step. Uses internal `last_order` tracking to count
consecutive steps at the same order.

Arguments:
- `stald`: StabilityLimitDetectionState
- `order`: the order used for the current step (NOT the newly chosen order)

Returns true if a stability violation is detected and order should be reduced by 1.
"""
function stald_check!(stald::StabilityLimitDetectionState, order::Int)
    if !stald.enabled
        return false
    end

    # Only check for orders >= 3 (1-2 are A-stable)
    if order < 3
        stald.last_order = order
        return false
    end

    # Track consecutive steps at the same order
    if order != stald.last_order
        stald.nscon = 0
        stald.last_order = order
        return false
    else
        stald.nscon += 1
    end

    # Need enough data: nscon >= order + 5
    if stald.nscon < order + 5
        return false
    end

    # Run the stability limit detection
    ldflag = _stald_detect(stald, order)

    if ldflag > 3
        # Stability violation detected
        stald.nor += 1
        stald.nscon = 0  # Reset counter after reduction
        return true
    end

    return false
end

"""
    stald_reset!(stald)

Reset STALD state (e.g., after u_modified or reinit).
"""
function stald_reset!(stald::StabilityLimitDetectionState{T}) where {T}
    z = zero(T)
    fill!(stald.ssdat, z)
    stald.nscon = 0
    stald.last_order = 0
    fill!(stald.smax, z)
    fill!(stald.ssmax, z)
    fill!(stald.rav, z)
    fill!(stald.vrat, z)
    fill!(stald.rat, z)
    fill!(stald.qc, z)
    fill!(stald.qco, z)
    fill!(stald.drr, z)
    fill!(stald.rrc, z)
    fill!(stald.sigsq, z)
    return nothing
end

"""
    _stald_detect(stald, q) -> Int

Core stability limit detection algorithm (port of CVODE's cvSLdet).

Returns:
- Negative values: algorithm could not determine (insufficient data quality)
- 1, 2, 3: characteristic root found, but rr <= rrcut (stable)
- 4, 5, 6: stability violation detected (rr > rrcut)
"""
function _stald_detect(stald::StabilityLimitDetectionState{T}, q::Int) where {T}
    ssdat = stald.ssdat
    (; rrcut, vrrtol, vrrt2, sqtol, rrtol, tiny) = stald
    (; smax, ssmax, rav, vrat, rat, qc, qco, drr, rrc, sigsq) = stald

    # Phase 1: Compute statistics from stored data
    # For each of the three derivative levels k=1,2,3
    z = zero(T)
    fill!(smax, z)
    fill!(ssmax, z)
    fill!(rav, z)
    fill!(vrat, z)
    fill!(rat, z)
    fill!(qc, z)

    for k in 1:3
        smink = ssdat[1, k]
        smaxk = z
        for i in 1:5
            smink = min(smink, ssdat[i, k])
            smaxk = max(smaxk, ssdat[i, k])
        end

        # Data too spread out
        if smink < tiny * smaxk
            return -1
        end

        smax[k] = smaxk
        ssmax[k] = smaxk * smaxk

        # Compute ratios of consecutive data points and their variance
        sumrat = z
        sumrsq = z
        for i in 1:4
            rat[i, k] = ssdat[i, k] / ssdat[i + 1, k]
            sumrat += rat[i, k]
            sumrsq += rat[i, k] * rat[i, k]
        end
        rav[k] = T(0.25) * sumrat                           # average ratio
        vrat[k] = abs(T(0.25) * sumrsq - rav[k] * rav[k])   # variance

        # Quartic polynomial coefficients
        qc[5, k] = ssdat[1, k] * ssdat[3, k] - ssdat[2, k] * ssdat[2, k]
        qc[4, k] = ssdat[2, k] * ssdat[3, k] - ssdat[1, k] * ssdat[4, k]
        qc[3, k] = z
        qc[2, k] = ssdat[2, k] * ssdat[5, k] - ssdat[3, k] * ssdat[4, k]
        qc[1, k] = ssdat[4, k] * ssdat[4, k] - ssdat[3, k] * ssdat[5, k]
    end

    # Phase 2: Determine rr (characteristic root magnitude)
    vmin = min(vrat[1], vrat[2], vrat[3])
    vmax = max(vrat[1], vrat[2], vrat[3])
    rr = z
    kflag = 0

    if vmin < vrrtol * vrrtol
        # Low variance case: ratios are consistent
        if vmax > vrrt2 * vrrt2
            return -2  # Some ratios are inconsistent
        end

        # Average the three ratio estimates
        rr = (rav[1] + rav[2] + rav[3]) / T(3)
        drrmax = z
        for k in 1:3
            adrr = abs(rav[k] - rr)
            drrmax = max(drrmax, adrr)
        end
        if drrmax > vrrt2
            return -3
        end
        kflag = 1  # Found root via normal matrix case
    else
        # Higher variance: use quartic method
        # Gaussian elimination on quartic coefficients
        copyto!(qco, qc)

        # Simple Gaussian elimination (matching CVODE logic)
        if abs(qco[4, 1]) > tiny * smax[1] * smax[1]
            r1 = qco[4, 2] / qco[4, 1]
            r2 = qco[4, 3] / qco[4, 1]
            for j in 5:-1:1
                if j != 4
                    qco[j, 2] = qco[j, 2] - r1 * qco[j, 1]
                    qco[j, 3] = qco[j, 3] - r2 * qco[j, 1]
                end
            end
        else
            return -4
        end

        if abs(qco[5, 2]) > tiny * smax[2] * smax[2]
            r1 = qco[5, 3] / qco[5, 2]
            for j in 5:-1:1
                if j != 4 && j != 5
                    qco[j, 3] = qco[j, 3] - r1 * qco[j, 2]
                end
            end
        else
            # Try alternate elimination
            if abs(qco[5, 3]) > tiny * smax[3] * smax[3]
                return -4
            end
        end

        # Solve for rr
        if abs(qco[4, 3]) < tiny
            return -4
        end
        rr = -qco[5, 3] / qco[4, 3]

        if rr < tiny || rr > T(100)
            return -5
        end

        # Verify rr satisfies all three quartics
        sqmax = z
        for k in 1:3
            qkr = qc[5, k] + rr * (qc[4, k] + rr * rr * (qc[2, k] + rr * qc[1, k]))
            saqk = abs(qkr) / ssmax[k]
            sqmax = max(sqmax, saqk)
        end

        if sqmax < sqtol
            kflag = 2  # Found root via quartic
        else
            # Newton corrections to improve rr
            kmin = 0
            for it in 1:3
                fill!(drr, z)
                fill!(rrc, z)

                for k in 1:3
                    qkr = qc[5, k] +
                        rr * (qc[4, k] + rr * rr * (qc[2, k] + rr * qc[1, k]))
                    qp = qc[4, k] +
                        rr * rr * (T(3) * qc[2, k] + rr * T(4) * qc[1, k])
                    if abs(qp) > tiny * ssmax[k]
                        drr[k] = -qkr / qp
                    end
                    rrc[k] = rr + drr[k]
                end

                # Pick correction giving smallest residual
                sqmin = T(Inf)
                for k in 1:3
                    if rrc[k] > 0
                        qkr_k = qc[5, k] + rrc[k] *
                            (
                            qc[4, k] +
                                rrc[k] * rrc[k] *
                                (qc[2, k] + rrc[k] * qc[1, k])
                        )
                        saqk = abs(qkr_k) / ssmax[k]
                        if saqk < sqmin
                            sqmin = saqk
                            kmin = k
                        end
                    end
                end

                if kmin == 0
                    return -6
                end
                rr = rrc[kmin]

                if sqmin < sqtol
                    kflag = 3  # Found root via Newton-corrected quartic
                    break
                end
            end

            if kflag == 0
                return -6  # Newton didn't converge
            end
        end
    end

    # Phase 3: Compute sigsq and cross-check rr
    fill!(sigsq, z)
    for k in 1:3
        rsa = ssdat[1, k]
        rsb = ssdat[2, k] * rr
        rsc = ssdat[3, k] * rr * rr
        rsd = ssdat[4, k] * rr * rr * rr
        rd1a = rsa - rsb
        rd1b = rsb - rsc
        rd1c = rsc - rsd
        rd2a = rd1a - rd1b
        rd2b = rd1b - rd1c
        rd3a = rd2a - rd2b

        if abs(rd1b) < tiny * smax[k]
            return -7
        end

        cest1 = -rd3a / rd1b
        if cest1 < tiny || cest1 > T(4)
            return -7
        end
        corr1 = (rd2b / cest1) / (rr * rr)
        sigsq[k] = ssdat[3, k] + corr1
    end

    if sigsq[2] < tiny
        return -8
    end

    # Phase 4: Cross-check rr from sigsq ratios
    ratp = sigsq[3] / sigsq[2]
    ratm = sigsq[1] / sigsq[2]
    qfac1 = T(0.25) * (q * q - 1)
    qfac2 = T(2) / (q - 1)
    bb = ratp * ratm - one(T) - qfac1 * ratp
    tem = one(T) - qfac2 * bb

    if abs(tem) < tiny
        return -8
    end

    rrb = one(T) / tem

    if abs(rrb - rr) > rrtol
        return -9
    end

    # Phase 5: Final stability decision
    if rr > rrcut
        # Stability violation detected!
        return kflag + 3  # Returns 4, 5, or 6
    end

    return kflag  # Returns 1, 2, or 3 (stable)
end
