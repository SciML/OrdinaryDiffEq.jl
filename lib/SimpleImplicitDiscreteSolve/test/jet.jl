import SimpleImplicitDiscreteSolve
using JET

@testset "JET Tests" begin
    test_package(
        SimpleImplicitDiscreteSolve, target_defined_modules = true, mode = :typo
    )
end
