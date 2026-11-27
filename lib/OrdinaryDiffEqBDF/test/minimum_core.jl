# Standalone release check: julia --startup-file=no minimum_core.jl REV [--git-deps]
# Before the Core and SDIRK releases, --git-deps takes both from the BDF revision.
using Pkg, Test

const repo = "https://github.com/SciML/OrdinaryDiffEq.jl.git"
gitpkg(name) = PackageSpec(; url = repo, rev = ARGS[1], subdir = "lib/$name")
minimum_pkg(name, version) = "--git-deps" in ARGS ? gitpkg(name) : PackageSpec(; name, version)
const bdf = gitpkg("OrdinaryDiffEqBDF")
const sdirk = minimum_pkg("OrdinaryDiffEqSDIRK", v"2.10.0")

Pkg.activate(; temp = true)
try
    Pkg.add([bdf, sdirk, PackageSpec(; name = "OrdinaryDiffEqCore", version = v"4.19.1")])
    error("BDF resolved with pre-hook Core")
catch err
    err isa Pkg.Resolve.ResolverError || rethrow()
    message = sprint(showerror, err)
    @test occursin("OrdinaryDiffEqCore", message)
    @test occursin("4.19.1", message)
    @test occursin("4.21", message)
    println(message)
end

Pkg.activate(; temp = true)
Pkg.add([bdf, sdirk, minimum_pkg("OrdinaryDiffEqCore", v"4.21.0")])
@test !any(info -> info.is_tracking_path, values(Pkg.dependencies()))
using OrdinaryDiffEqBDF
for prob in (
        ODEProblem((u, p, t) -> -u, 1.0, (0.0, 0.1)),
        ODEProblem((du, u, p, t) -> (du[1] = -u[1]), [1.0], (0.0, 0.1)),
    )
    dt = init(prob, FBDF(); abstol = 1.0e-8, reltol = 1.0e-6).dt
    @show dt
    @test dt ≈ sqrt(0.01 * (1.0e-8 + 1.0e-6)) rtol = 100eps()
end
