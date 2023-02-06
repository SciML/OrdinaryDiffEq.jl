function explicit_rk_docstring(description::String)
    start_docstring = """
       ```julia
       $FUNCTIONNAME(; stage_limiter! = OrdinaryDiffEq.trivial_limiter!,
       step_limiter! = OrdinaryDiffEq.trivial_limiter!,
       thread = OrdinaryDiffEq.False())
       ```

       Explicit Runge-Kutta Method.
       """
    end_docstring = """

        ### Keyword Arguments

         - `stage_limiter!`: function of the form `limiter!(u, integrator, p, t)`
         - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
         - `thread`: determines whether internal broadcasting on
            appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`,
            default) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when
            Julia is started with multiple threads.
        """
    start_docstring * description * end_docstring
end

@doc explicit_rk_docstring("The second order Heun's method. Uses embedded Euler method for adaptivity.")
struct Heun{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function Heun(; stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
              thread = False())
    Heun{typeof(stage_limiter!), typeof(step_limiter!), typeof(thread)}(stage_limiter!,
                                                                        step_limiter!,
                                                                        thread)
end

# for backwards compatibility
function Heun(stage_limiter!, step_limiter! = trivial_limiter!)
    Heun{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
                                                               step_limiter!,
                                                               False())
end

function Base.show(io::IO, alg::Heun)
    print(io, "Heun(stage_limiter! = ", alg.stage_limiter!,
          ", step_limiter! = ", alg.step_limiter!,
          ", thread = ", alg.thread, ")")
end
