alg_order(alg::RandomEM) = 1 // 2
alg_order(alg::RandomHeun) = 1 // 2
alg_order(alg::RandomTamedEM) = 1 // 2
alg_order(alg::RandomTaylor15) = 3 // 2
alg_order(alg::BAOAB) = 1 // 1

alg_compatible(prob::SciMLBase.AbstractSDEProblem, alg::BAOAB) = is_diagonal_noise(prob)

function alg_compatible(prob::SciMLBase.AbstractRODEProblem, alg::RandomTaylor15)
    return prob.noise isa DiffEqNoiseProcess.NoiseGrid && eltype(prob.noise.W) <: Number
end
