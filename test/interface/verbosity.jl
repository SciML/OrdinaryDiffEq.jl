using ODEProblemLibrary: prob_ode_vanderpol_stiff
using OrdinaryDiffEq, SciMLBase
import SciMLBase: Verbosity, ODEPerformanceVerbosity, ODEVerbosity


# Stiff Switching

verb = ODEVerbosity(performance = ODEPerformanceVerbosity(alg_switch = Verbosity.Info()))
solve(prob_ode_vanderpol_stiff, AutoTsit5(Rodas5()), verbose = verb)