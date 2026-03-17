# convenience macro
macro wrap_h(signature)
    Meta.isexpr(signature, :call) ||
        throw(ArgumentError("signature has to be a function call expression"))

    name = signature.args[1]
    args = signature.args[2:end]
    args_wo_h = [arg for arg in args if arg !== :h]

    return quote
        if f.$name === nothing
            nothing
        else
            if isinplace(f)
                let _f = f.$name, h = h
                    ($(args_wo_h...),) -> _f($(args...))
                end
            else
                let _f = f.$name, h = h
                    ($(args_wo_h[2:end]...),) -> _f($(args[2:end]...))
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

function ODEFunctionWrapper(f::SciMLBase.AbstractDDEFunction, h)
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
