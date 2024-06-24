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