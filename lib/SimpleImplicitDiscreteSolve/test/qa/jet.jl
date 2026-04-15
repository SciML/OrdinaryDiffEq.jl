import SimpleImplicitDiscreteSolve
using JET

@testset "JET Tests" begin
    test_package(
        SimpleImplicitDiscreteSolve, target_modules = (SimpleImplicitDiscreteSolve,), mode = :typo
    )
end
