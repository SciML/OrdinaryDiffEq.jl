using Pkg
Pkg.add("JET")

import OrdinaryDiffEqAdamsBashforthMoulton
using JET

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqAdamsBashforthMoulton, target_modules = (OrdinaryDiffEqAdamsBashforthMoulton,), mode = :typo
    )
end
