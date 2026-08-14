using DiffEqBase, SciMLBase, Test

if isdefined(SciMLBase, :report_integrator_failure)
    @test any(m -> m.module === DiffEqBase, methods(SciMLBase.report_integrator_failure))
else
    @test !isdefined(SciMLBase, :report_integrator_failure)
end
