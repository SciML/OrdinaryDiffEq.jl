# Linear Methods

for Alg in [
    :MagnusMidpoint,
    :MagnusLeapfrog,
    :LieEuler,
    :MagnusGauss4,
    :MagnusNC6,
    :MagnusGL6,
    :MagnusGL8,
    :MagnusNC8,
    :MagnusGL4,
    :RKMK2,
    :RKMK4,
    :LieRK4,
    :CG2,
    :CG3,
    :CG4a
]
    @eval struct $Alg <: OrdinaryDiffEqLinearExponentialAlgorithm
        krylov::Bool
        m::Int
        iop::Int
    end
    @eval $Alg(; krylov = false, m = 30, iop = 0) = $Alg(krylov, m, iop)
end

struct MagnusAdapt4 <: OrdinaryDiffEqAdaptiveAlgorithm end

struct CayleyEuler <: OrdinaryDiffEqAlgorithm end

struct LinearExponential <:
       OrdinaryDiffEqExponentialAlgorithm{1, false, Val{:forward}, Val{true}, nothing}
    krylov::Symbol
    m::Int
    iop::Int
end
LinearExponential(; krylov = :off, m = 10, iop = 0) = LinearExponential(krylov, m, iop)