using Pkg
Pkg.add("JET")

import OrdinaryDiffEqIMEXMultistep
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqIMEXMultistep, mode = :typo
    )
end
