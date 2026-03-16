# Backward-compat: LinearInterpolationData is used by StochasticDelayDiffEq.jl.
# SDE solutions now use ODE's isempty(k) linear interpolation fallback, but
# SDDE creates LinearInterpolationData directly as a custom solution interp.

struct LinearInterpolationData{uType, tType} <: DiffEqBase.AbstractDiffEqInterpolation
    timeseries::uType
    ts::tType
end

DiffEqBase.interp_summary(::LinearInterpolationData) = "1st order linear"

function (interp::LinearInterpolationData)(
        tvals, idxs, deriv, p, continuity::Symbol = :left
    )
    return _sde_linear_interpolation(tvals, interp, idxs, deriv, continuity)
end

function (interp::LinearInterpolationData)(
        val, tvals, idxs, deriv, p, continuity::Symbol = :left
    )
    return _sde_linear_interpolation!(val, tvals, interp, idxs, deriv, continuity)
end

# --- out-of-place ---

function _sde_linear_interpolation(
        tvals, id, idxs, deriv, continuity::Symbol
    )
    (; ts, timeseries) = id
    @inbounds tdir = sign(ts[end] - ts[1])
    idx = sortperm(tvals, rev = tdir < 0)

    i₋ = 1
    i₊ = 2
    vals = map(idx) do j
        @inbounds begin
            t = tvals[j]
            if continuity === :left
                i₊ = min(
                    lastindex(ts), OrdinaryDiffEqCore._searchsortedfirst(
                        ts, t, i₊, tdir > 0
                    )
                )
                i₋ = i₊ > 1 ? i₊ - 1 : i₊
            else
                i₋ = max(
                    1, OrdinaryDiffEqCore._searchsortedlast(ts, t, i₋, tdir > 0)
                )
                i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
            end
            dt = ts[i₊] - ts[i₋]
            Θ = iszero(dt) ? oneunit(t) / oneunit(dt) : (t - ts[i₋]) / dt
            _sde_linear_interpolant(Θ, dt, timeseries[i₋], timeseries[i₊], idxs, deriv)
        end
    end
    tdir < 0 && permute!(vals, idx)
    return DiffEqBase.DiffEqArray(vals, tvals)
end

function _sde_linear_interpolation(
        tval::Number, id, idxs, deriv, continuity::Symbol
    )
    (; ts, timeseries) = id
    @inbounds tdir = sign(ts[end] - ts[1])
    if continuity === :left
        i₊ = min(
            lastindex(ts), OrdinaryDiffEqCore._searchsortedfirst(
                ts, tval, 2, tdir > 0
            )
        )
        i₋ = i₊ > 1 ? i₊ - 1 : i₊
    else
        i₋ = max(1, OrdinaryDiffEqCore._searchsortedlast(ts, tval, 1, tdir > 0))
        i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
    end
    @inbounds begin
        dt = ts[i₊] - ts[i₋]
        Θ = iszero(dt) ? oneunit(tval) / oneunit(dt) : (tval - ts[i₋]) / dt
        return _sde_linear_interpolant(
            Θ, dt, timeseries[i₋], timeseries[i₊], idxs, deriv
        )
    end
end

# --- in-place ---

function _sde_linear_interpolation!(
        out, tval::Number, id, idxs, deriv, continuity::Symbol
    )
    (; ts, timeseries) = id
    @inbounds tdir = sign(ts[end] - ts[1])
    if continuity === :left
        i₊ = min(
            lastindex(ts), OrdinaryDiffEqCore._searchsortedfirst(
                ts, tval, 2, tdir > 0
            )
        )
        i₋ = i₊ > 1 ? i₊ - 1 : i₊
    else
        i₋ = max(1, OrdinaryDiffEqCore._searchsortedlast(ts, tval, 1, tdir > 0))
        i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
    end
    @inbounds begin
        dt = ts[i₊] - ts[i₋]
        Θ = iszero(dt) ? oneunit(tval) / oneunit(dt) : (tval - ts[i₋]) / dt
        _sde_linear_interpolant!(
            out, Θ, dt, timeseries[i₋], timeseries[i₊], idxs, deriv
        )
    end
    return out
end

# --- low-level linear interpolants ---

@muladd @inline function _sde_linear_interpolant(
        Θ, dt, u0, u1, idxs::Nothing, ::Type{Val{0}}
    )
    @.. (1 - Θ) * u0 + Θ * u1
end

@muladd @inline function _sde_linear_interpolant(
        Θ, dt, u0, u1, idxs, ::Type{Val{0}}
    )
    @.. (1 - Θ) * u0[idxs] + Θ * u1[idxs]
end

@inline function _sde_linear_interpolant(
        Θ, dt, u0, u1, idxs::Nothing, ::Type{Val{1}}
    )
    @.. (u1 - u0) / dt
end

@inline function _sde_linear_interpolant(Θ, dt, u0, u1, idxs, ::Type{Val{1}})
    @.. (u1[idxs] - u0[idxs]) / dt
end

@muladd @inline function _sde_linear_interpolant!(
        out, Θ, dt, u0, u1, idxs, ::Type{Val{0}}
    )
    Θm1 = (1 - Θ)
    if idxs === nothing
        @.. out = Θm1 * u0 + Θ * u1
    else
        @views @.. out = Θm1 * u0[idxs] + Θ * u1[idxs]
    end
end

@inline function _sde_linear_interpolant!(
        out, Θ, dt, u0, u1, idxs, ::Type{Val{1}}
    )
    if idxs === nothing
        @.. out = (u1 - u0) / dt
    else
        @views @.. out = (u1[idxs] - u0[idxs]) / dt
    end
end
