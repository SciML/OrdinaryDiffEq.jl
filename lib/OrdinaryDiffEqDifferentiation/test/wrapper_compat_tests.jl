using OrdinaryDiffEqDifferentiation, Test

const SciMLBaseOwner = OrdinaryDiffEqDifferentiation.SciMLBase

@test OrdinaryDiffEqDifferentiation.TimeDerivativeWrapper ===
    SciMLBaseOwner.TimeDerivativeWrapper
@test OrdinaryDiffEqDifferentiation.TimeGradientWrapper ===
    SciMLBaseOwner.TimeGradientWrapper
@test OrdinaryDiffEqDifferentiation.UDerivativeWrapper ===
    SciMLBaseOwner.UDerivativeWrapper
@test OrdinaryDiffEqDifferentiation.UJacobianWrapper ===
    SciMLBaseOwner.UJacobianWrapper
