# convenience macro
macro wrap_h(signature)
    Meta.isexpr(signature, :call) ||
        throw(ArgumentError("signature has to be a function call expression"))

    name = signature.args[1]
    args = signature.args[2:end]
    args_wo_h = [arg for arg in args if arg !== :h]

    # NOTE: we build the lambda parameter list explicitly as an `Expr(:tuple, …)`
    # rather than writing `($(args_wo_h...))` directly. Writing
    # `($(args_wo_h...)) -> body` in an interpolated quote with N > 1 elements
    # parses as `Expr(:->, :a, :b, …, body)` (which Julia then interprets as a
    # 1-argument destructuring lambda), not as a multi-argument lambda. The
    # explicit `Expr(:tuple, …)` avoids that pitfall and is Runic-stable.
    ip_args = Expr(:tuple, args_wo_h...)
    ip_call = Expr(:call, :_f, args...)
    oop_args = Expr(:tuple, args_wo_h[2:end]...)
    oop_call = Expr(:call, :_f, args[2:end]...)

    return quote
        if f.$name === nothing
            nothing
        else
            if isinplace(f)
                let _f = f.$name, h = h
                    $(Expr(:->, ip_args, ip_call))
                end
            else
                let _f = f.$name, h = h
                    $(Expr(:->, oop_args, oop_call))
                end
            end
        end
    end |> esc
end

struct ODEFunctionWrapper{iip, F, H, TMM, Ta, Tt, TJ, JP, SP, TW, TWt, TPJ, S, TCV, ID} <:
    SciMLBase.AbstractODEFunction{iip}
    f::F
    h::H
    mass_matrix::TMM
    analytic::Ta
    tgrad::Tt
    jac::TJ
    jac_prototype::JP
    sparsity::SP
    Wfact::TW
    Wfact_t::TWt
    paramjac::TPJ
    sys::S
    colorvec::TCV
    initialization_data::ID
end

function ODEFunctionWrapper(f::Union{SciMLBase.AbstractDDEFunction, SciMLBase.AbstractSDDEFunction}, h)
    # wrap functions
    jac = @wrap_h jac(J, u, h, p, t)
    Wfact = @wrap_h Wfact(W, u, h, p, dtgamma, t)
    Wfact_t = @wrap_h Wfact_t(W, u, h, p, dtgamma, t)

    return ODEFunctionWrapper{
        isinplace(f), typeof(f.f), typeof(h), typeof(f.mass_matrix),
        typeof(f.analytic), typeof(f.tgrad), typeof(jac),
        typeof(f.jac_prototype), typeof(f.sparsity),
        typeof(Wfact), typeof(Wfact_t),
        typeof(f.paramjac), typeof(f.sys), typeof(f.colorvec),
        typeof(f.initialization_data),
    }(
        f.f, h,
        f.mass_matrix,
        f.analytic,
        f.tgrad, jac,
        f.jac_prototype,
        f.sparsity,
        Wfact,
        Wfact_t,
        f.paramjac,
        f.sys,
        f.colorvec,
        f.initialization_data
    )
end

(f::ODEFunctionWrapper{true})(du, u, p, t) = f.f(du, u, f.h, p, t)
(f::ODEFunctionWrapper{false})(u, p, t) = f.f(u, f.h, p, t)

# Wrapper for the diffusion function g with history injection
struct DiffusionWrapper{iip, G, H}
    g::G
    h::H
end

(gw::DiffusionWrapper{true})(du, u, p, t) = gw.g(du, u, gw.h, p, t)
(gw::DiffusionWrapper{false})(u, p, t) = gw.g(u, gw.h, p, t)

"""
    SDEFunctionWrapper

Wraps an SDDEFunction (with history-dependent drift `f` and diffusion `g`)
to produce an SDE-like function interface `f(du, u, p, t)` and `f.g(du, u, p, t)`
that SDE algorithms can call. The history function `h` is captured in closures.
"""
struct SDEFunctionWrapper{iip, F, G, H, TMM, Ta, Tt, TJ, JP, SP, TW, TWt, TPJ, GG, S, TCV, ID} <:
    SciMLBase.AbstractSDEFunction{iip}
    f::F
    g::G
    h::H
    mass_matrix::TMM
    analytic::Ta
    tgrad::Tt
    jac::TJ
    jac_prototype::JP
    sparsity::SP
    Wfact::TW
    Wfact_t::TWt
    paramjac::TPJ
    ggprime::GG
    sys::S
    colorvec::TCV
    initialization_data::ID
end

function SDEFunctionWrapper(f::SciMLBase.AbstractSDDEFunction, h)
    # wrap jacobian-related functions
    jac = @wrap_h jac(J, u, h, p, t)
    Wfact = @wrap_h Wfact(W, u, h, p, dtgamma, t)
    Wfact_t = @wrap_h Wfact_t(W, u, h, p, dtgamma, t)

    # wrap diffusion with history
    g_wrapped = DiffusionWrapper{isinplace(f), typeof(f.g), typeof(h)}(f.g, h)

    return SDEFunctionWrapper{
        isinplace(f), typeof(f.f), typeof(g_wrapped), typeof(h),
        typeof(f.mass_matrix),
        typeof(f.analytic), typeof(f.tgrad), typeof(jac),
        typeof(f.jac_prototype), typeof(f.sparsity),
        typeof(Wfact), typeof(Wfact_t),
        typeof(f.paramjac),
        typeof(f.ggprime),
        typeof(f.sys), typeof(f.colorvec),
        typeof(f.initialization_data),
    }(
        f.f, g_wrapped, h,
        f.mass_matrix,
        f.analytic,
        f.tgrad, jac,
        f.jac_prototype,
        f.sparsity,
        Wfact,
        Wfact_t,
        f.paramjac,
        f.ggprime,
        f.sys,
        f.colorvec,
        f.initialization_data
    )
end

(f::SDEFunctionWrapper{true})(du, u, p, t) = f.f(du, u, f.h, p, t)
(f::SDEFunctionWrapper{false})(u, p, t) = f.f(u, f.h, p, t)
