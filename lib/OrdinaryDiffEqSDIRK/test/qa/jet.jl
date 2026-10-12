import OrdinaryDiffEqSDIRK
using JET

@testset "JET Tests" begin
    if VERSION >= v"1.12"
        # JET 0.12 reports the stage variables of the unrolled out-of-place IMEX step, which are
        # assigned under `s >= n` guards, as possibly undefined.
        known(report) = report isa JET.UndefVarErrorReport && any(
            frame -> frame.linfo isa Core.MethodInstance &&
                frame.linfo.def.name === :_perform_step_oop!,
            report.vst
        )
        result = JET.report_package(
            OrdinaryDiffEqSDIRK; target_modules = (OrdinaryDiffEqSDIRK,), mode = :typo
        )
        @test isempty(filter(!known, JET.get_reports(result)))
    else
        test_package(
            OrdinaryDiffEqSDIRK, target_modules = (OrdinaryDiffEqSDIRK,), mode = :typo
        )
    end
end
