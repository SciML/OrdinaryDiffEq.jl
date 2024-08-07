"""
Utility function to help generating consistent docstrings across the package.
"""
# TODO we should add a consistency check using the string $name(; $keyword_default) against the standard console output of name()
function generic_solver_docstring(description::String,
        name::String,
        solver_class::String,
        references::String,
        keyword_description::String,
        keyword_default::String)
    if !isempty(keyword_default)
        # Chunk string and remove empty lines
        kws = split(keyword_default, "\n")
        keywords_split = [kw for kw in kws if !isempty(rstrip(kw))]

        # Indent the keywords properly
        indentation = repeat(" ", length(name) + 3)
        # We do not indent the first kw and no newline for the last one
        if length(keywords_split) > 1
            keywords_split[1] = keywords_split[1] * "\n"
            for i in 2:(length(keywords_split) - 1)
                keywords_split[i] = indentation * keywords_split[i] * "\n"
            end
            keywords_split[end] = indentation * keywords_split[end]
        end
        # Flatten into string
        keyword_default = join(keywords_split)
        # Remove trailing comma
        keyword_default = strip(keyword_default, [','])
    end
    start_docstring = !isempty(keyword_default) ? """
        ```julia
        $name(; $keyword_default)
        ```

        $solver_class
        """ : """
        ```julia
        $name()
        ```
        $solver_class
        """

    keyword_docstring = """

        ### Keyword Arguments

        $keyword_description
        """

    return start_docstring * description * keyword_docstring *
           "## References\n" * references
end

function explicit_rk_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description::String = "",
        extra_keyword_default::String = "")
    keyword_default = """
    stage_limiter! = OrdinaryDiffEq.trivial_limiter!,
    step_limiter! = OrdinaryDiffEq.trivial_limiter!,
    thread = OrdinaryDiffEq.False(),
    """ * extra_keyword_default

    keyword_default_description = """
    - `stage_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    - `thread`: determines whether internal broadcasting on appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads.
    """ * extra_keyword_description

    generic_solver_docstring(
        description, name, "Explicit Runge-Kutta Method. ", references,
        keyword_default_description, keyword_default
    )
end
function differentiation_rk_docstring(description::String,
        name::String,
        solver_class::String;
        references::String = "",
        extra_keyword_description::String = "",
        extra_keyword_default::String = "")
    keyword_default = """
    chunk_size = Val{0}(),
    autodiff = true,
    standardtag = Val{true}(),
    concrete_jac = nothing,
    diff_type = Val{:forward},
    linsolve = nothing,
    """ * extra_keyword_default

    keyword_default_description = """
    - `chunk_size`: The chunk size used with ForwardDiff.jl. Defaults to `Val{0}()`
        and thus uses the internal ForwardDiff.jl algorithm for the choice.
    - `autodiff`: Specifies whether to use automatic differentiation via
        [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) or finite
        differencing via [FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl).
        Defaults to `Val{true}()` for automatic differentiation.
    - `standardtag`: Specifies whether to use package-specific tags instead of the
        ForwardDiff default function-specific tags. For more information, see
        [this blog post](https://www.stochasticlifestyle.com/improved-forwarddiff-jl-stacktraces-with-package-tags/).
        Defaults to `Val{true}()`.
    - `concrete_jac`: Specifies whether a Jacobian should be constructed. Defaults to
        `nothing`, which means it will be chosen true/false depending on circumstances
        of the solver, such as whether a Krylov subspace method is used for `linsolve`.
    - `diff_type`: The type of differentiation used in FiniteDiff.jl if `autodiff=false`.
        Defaults to `Val{:forward}`, with alternatives of `Val{:central}` and
        `Val{:complex}`.
    - `linsolve`: Any [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) compatible linear solver.
      For example, to use [KLU.jl](https://github.com/JuliaSparse/KLU.jl), specify
      `$name(linsolve = KLUFactorization()`).
       When `nothing` is passed, uses `DefaultLinearSolver`.

    """ * extra_keyword_description

    generic_solver_docstring(
        description, name, solver_class, references,
        keyword_default_description, keyword_default
    )
end
