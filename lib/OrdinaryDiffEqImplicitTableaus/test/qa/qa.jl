using OrdinaryDiffEqImplicitTableaus
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(OrdinaryDiffEqImplicitTableaus)
end

@testset "JET" begin
    JET.test_package(OrdinaryDiffEqImplicitTableaus; target_defined_modules = true)
end
