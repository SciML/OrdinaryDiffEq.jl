import OrdinaryDiffEqSDIRK
using JET

@testset "JET Tests" begin
    if VERSION >= v"1.12"
        # JET 0.12 reports the stage variables of the unrolled out-of-place IMEX step, which are
        # assigned under `s >= n` guards, as possibly undefined. Only those locals are skipped.
        function known(report)
            report isa JET.UndefVarErrorReport && report.maybeundef || return false
            report.var isa Symbol && occursin(r"^([zk]\d+|tmp_est)$", String(report.var)) ||
                return false
            m = report.vst[end].linfo.def
            return m isa Method && m.module === OrdinaryDiffEqSDIRK &&
                m.name === :_perform_step_oop! &&
                basename(String(m.file)) == "generic_imex_perform_step.jl"
        end
        result = JET.report_package(
            OrdinaryDiffEqSDIRK; target_modules = (OrdinaryDiffEqSDIRK,), mode = :typo
        )
        @test isempty(filter(!known, JET.get_reports(result)))
        # A real typo in a `_perform_step_oop!` method must still be reported.
        Core.eval(OrdinaryDiffEqSDIRK, :(_perform_step_oop!(::Val{:jet_control}) = reviewer_misspelled_stage))
        control = JET.get_reports(
            JET.report_call(OrdinaryDiffEqSDIRK._perform_step_oop!, (Val{:jet_control},); mode = :typo)
        )
        @test !isempty(control) && !any(known, control)
    else
        test_package(
            OrdinaryDiffEqSDIRK, target_modules = (OrdinaryDiffEqSDIRK,), mode = :typo
        )
    end
end
