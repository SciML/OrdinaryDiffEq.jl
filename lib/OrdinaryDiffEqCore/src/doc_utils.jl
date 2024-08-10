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
        if length(keyword_default) > 1
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
        """ :
                      """
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
    """ * extra_keyword_default

    keyword_default_description = """
    - `stage_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """ * extra_keyword_description

    generic_solver_docstring(
        description, name, "Explicit Runge-Kutta Method. ", references,
        keyword_default_description, keyword_default
    )
end
