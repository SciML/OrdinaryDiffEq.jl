function isdtchangeable(alg::Union{LawsonEuler, NorsettEuler, ETDRK2, ETDRK3, ETDRK4, HochOst4, ETD2})
    false
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

function DiffEqBase.prepare_alg(
        alg::ETD2,
        u0::AbstractArray,
        p, prob)
    alg
end

fsal_typeof(alg::ETD2, rate_prototype) = ETD2Fsal{typeof(rate_prototype)}
function fsal_typeof(alg::CompositeAlgorithm, rate_prototype)
    fsal = map(x -> fsal_typeof(x, rate_prototype), alg.algs)
    @assert length(unique(fsal))==1 "`fsal_typeof` must be consistent"
    return fsal[1]
end

ismultistep(alg::ETD2) = true