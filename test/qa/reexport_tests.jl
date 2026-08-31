using SciMLBase, Test

@testset "Re-exported SciMLBase API follows real OrdinaryDiffEq usage" begin
    repo_dir = joinpath(@__DIR__, "..", "..")

    documented = Symbol[]
    for line in eachline(joinpath(repo_dir, "docs", "src", "api", "reexports.md"))
        occursin(r"^[A-Za-z_][A-Za-z0-9_!]*(?:, [A-Za-z_][A-Za-z0-9_!]*)*,?$", line) || continue
        append!(documented, Symbol.(filter(!isempty, strip.(split(line, ",")))))
    end
    @test allunique(documented)
    @test length(documented) > 50

    exported = setdiff(names(SciMLBase), [:SciMLBase])
    is_scimlbase_api(name) = @static if isdefined(Base, :ispublic)
        Base.ispublic(SciMLBase, name)
    else
        isdefined(SciMLBase, name)
    end
    @test all(is_scimlbase_api, documented)

    common_interface = read(
        joinpath(repo_dir, "docs", "src", "api", "common_interface.md"), String
    )
    referenced = Set(
        Symbol(m.captures[1]) for
            m in eachmatch(r"SciMLBase\.([A-Za-z_][A-Za-z0-9_!]*)", common_interface)
    )
    is_concrete_type(name) = begin
        value = getglobal(SciMLBase, name)
        value isa Type && !isabstracttype(Base.unwrap_unionall(value))
    end
    problem_api = filter(exported) do name
        text = String(name)
        (endswith(text, "Problem") || endswith(text, "Function")) &&
            is_concrete_type(name)
    end
    required_problem_api = intersect(problem_api, referenced)
    @test issetequal(intersect(documented, problem_api), required_problem_api)

    function reexported_names(path)
        source = read(path, String)
        m = match(r"@reexport using SciMLBase:((?:[^\n]*\n)(?:    [^\n]*\n)*)", source)
        m === nothing && return nothing
        return Symbol.(filter(!isempty, strip.(split(replace(m[1], r"\s+" => " "), ","))))
    end

    lib_dir = joinpath(repo_dir, "lib")
    selective = Dict{String, Vector{Symbol}}()
    for package in readdir(lib_dir)
        source = joinpath(lib_dir, package, "src", package * ".jl")
        isfile(source) || continue
        found = reexported_names(source)
        found === nothing || (selective[package] = found)
    end
    @test length(selective) == 9
    for (package, found) in selective
        @test (package, sort(found)) == (package, sort(documented))
    end
end
