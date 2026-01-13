using Pkg
Pkg.add("JET")

import OrdinaryDiffEqDefault
using JET

if isempty(VERSION.prerelease)
    @testset "JET Tests" begin
        test_package(
            OrdinaryDiffEqDefault, target_defined_modules = true, mode = :typo
        )
    end
end
