using ImplicitDiscreteSolve
using JET
using Test

@testset "JET Tests" begin
    test_package(
        ImplicitDiscreteSolve, target_modules = (ImplicitDiscreteSolve,), mode = :typo
    )
end
