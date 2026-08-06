function isdtchangeable(
        alg::Union{
            LawsonEuler, NorsettEuler, ETDRK2, ETDRK3, ETDRK4, HochOst4, ETD2,
        }
    )
    return false
end # due to caching

alg_order(alg::LawsonEuler) = 1
alg_order(alg::NorsettEuler) = 1
alg_order(alg::ETDRK2) = 2
alg_order(alg::ETDRK3) = 3
alg_order(alg::ETDRK4) = 4
alg_order(alg::HochOst4) = 4
alg_order(alg::Exp4) = 4
alg_order(alg::EPIRK4s3A) = 4
alg_order(alg::EPIRK4s3B) = 4
alg_order(alg::EPIRK5s3) = 5
alg_order(alg::EPIRK5P1) = 5
alg_order(alg::EPIRK5P2) = 5
alg_order(alg::EXPRB53s3) = 5
alg_order(alg::ETD2) = 2
alg_order(alg::Exprb32) = 3
alg_order(alg::Exprb43) = 4

alg_adaptive_order(alg::Exprb32) = 2
alg_adaptive_order(alg::Exprb43) = 4

function _expRK_has_concretization(A)
    return try
        has_concretization(A)
    catch err
        err isa UndefVarError || rethrow()
        true
    end
end

_expRK_requires_krylov(f) = false
function _expRK_requires_krylov(f::SplitFunction)
    A = f.f1.f
    return size(A) != () && !_expRK_has_concretization(A)
end

for Alg in (LawsonEuler, NorsettEuler, ETDRK2, ETDRK3, ETDRK4, HochOst4)
    @eval function DiffEqBase.prepare_alg(
            alg::$Alg,
            u0::AbstractArray,
            p, prob
        )
        if !alg.krylov && _expRK_requires_krylov(prob.f)
            return $Alg(
                krylov = true, m = alg.m, iop = alg.iop,
                autodiff = alg.autodiff, concrete_jac = alg.concrete_jac
            )
        end
        return alg
    end
end

function DiffEqBase.prepare_alg(
        alg::ETD2,
        u0::AbstractArray,
        p, prob
    )
    return alg
end

fsal_typeof(alg::ETD2, rate_prototype) = ETD2Fsal{typeof(rate_prototype)}

ismultistep(alg::ETD2) = true
