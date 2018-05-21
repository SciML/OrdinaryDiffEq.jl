using DiffEqBase
using OrdinaryDiffEq
using Base.Test

@test DiffEqBase.undefined_exports(OrdinaryDiffEq) == []
