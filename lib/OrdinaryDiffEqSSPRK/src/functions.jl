@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::SSPRK22,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end

@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::SSPRK33,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end

@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::SSPRK53_2N1,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end

@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::SSPRK53_2N2,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end

@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::SSPRK432,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end

@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::SSPRK932,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end

@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end

@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end

@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end

@eval @inline function DiffEqBase.get_tmp_cache(integrator, alg::OrdinaryDiffEqAlgorithm,
    cache::OrdinaryDiffEqConstantCache)
    nothing
end