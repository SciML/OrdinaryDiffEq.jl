using Pkg
Pkg.add("JET")

import SimpleImplicitDiscreteSolve
using JET

@testset "JET Tests" begin
    test_package(
        SimpleImplicitDiscreteSolve, mode = :typo
    )
end
