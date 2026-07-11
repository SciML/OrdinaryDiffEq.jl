using SciMLTesting, OrdinaryDiffEqDefault, Test
using Pkg: Pkg

# `DefaultSolverChoice` is an `EnumX.@enumx`-generated submodule that
# ExplicitImports cannot statically analyze. Its members (`Tsit5`, `Vern7`,
# `Rosenbrock23`, `Rodas5P`, `FBDF`, ...) share names with the solver types
# imported into `OrdinaryDiffEqDefault`, so the stale-import analysis
# mis-attributes those still-used imports as unused.
const ENUM_SUBMODULE = (OrdinaryDiffEqDefault.DefaultSolverChoice,)

function without_local_project_sources(f, pkg)
    project_file = joinpath(pkgdir(pkg), "Project.toml")
    manifest_file = joinpath(pkgdir(pkg), "Manifest.toml")
    project = Pkg.TOML.parsefile(project_file)
    if !haskey(project, "sources") && !isfile(manifest_file)
        return f()
    end
    original = read(project_file, String)
    manifest = isfile(manifest_file) ? read(manifest_file, String) : nothing
    delete!(project, "sources")
    open(project_file, "w") do io
        Pkg.TOML.print(io, project)
    end
    rm(manifest_file; force = true)
    try
        return f()
    finally
        open(project_file, "w") do io
            write(io, original)
        end
        if manifest !== nothing
            open(manifest_file, "w") do io
                write(io, manifest)
            end
        end
    end
end

without_local_project_sources(OrdinaryDiffEqDefault) do
    run_qa(
        OrdinaryDiffEqDefault;
        aqua_kwargs = (; piracies = false),
        explicit_imports = true,
        ei_kwargs = (;
            no_implicit_imports = (; allow_unanalyzable = ENUM_SUBMODULE),
            no_stale_explicit_imports = (;
                allow_unanalyzable = ENUM_SUBMODULE,
                # Used by `default_alg.jl` (solver constructors) but mis-flagged
                # stale because of the unanalyzable enum submodule above.
                ignore = (:Tsit5, :Vern7, :Rosenbrock23, :Rodas5P, :FBDF),
            ),
            all_explicit_imports_are_public = (;
                # Internal monorepo extension point; `AutoSwitchCache` is not user API.
                ignore = (:AutoSwitchCache,),
            ),
            all_qualified_accesses_are_public = (;
                # `lorenz`/`lorenz_oop` are OrdinaryDiffEqCore precompile-workload
                # test problems (defined in `precompilation_setup.jl`), deliberately
                # not part of its public extension surface.
                ignore = (:lorenz, :lorenz_oop),
            ),
        ),
        api_docs_kwargs = (;
            # Reexported upstream SciMLOperators names are documented at their owner.
            ignore = (:StaticWOperator, :has_concretization),
        ),
    )
end
