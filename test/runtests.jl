using Pkg
using SafeTestsets, Test
using SciMLTesting

const LONGER_TESTS = false
const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")

const GROUP = current_group()
const LIB_DIR = joinpath(dirname(@__DIR__), "lib")

# Folder-based v1.2 layout: each functional group's test files live in its own
# `test/<GroupKey>/` folder (matching the test_groups.toml key casing). The group
# bodies below `@safetestset include` the files from those folders, reproducing the
# corresponding `if GROUP == ...` branch of the previous hand-written runtests.jl
# verbatim (same files, same order, same @safetestset isolation, same version /
# is_APPVEYOR guards). They are passed to `run_tests` as 0-arg thunks so the version
# gates survive (folder-discovery mode cannot express them) and so the `@safetestset`
# macros resolve against the `using SafeTestsets` in this module's scope.

function interface_i()
    # Skip on Julia LTS (oneunit(Type{Any}) not defined) and pre-release (stalls)
    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/2979
    if VERSION >= v"1.11" && isempty(VERSION.prerelease)
        @time @safetestset "Null u0 Callbacks Tests" include("InterfaceI/null_u0_callbacks_test.jl")
    end
    @time @safetestset "Tstops Tests" include("InterfaceI/ode_tstops_tests.jl")
    @time @safetestset "Backwards Tests" include("InterfaceI/ode_backwards_test.jl")
    @time @safetestset "Initdt Tests" include("InterfaceI/ode_initdt_tests.jl")
    @time @safetestset "Linear Tests" include("InterfaceI/ode_twodimlinear_tests.jl")
    @time @safetestset "Inf Tests" include("InterfaceI/inf_handling.jl")
    @time @safetestset "saveat Tests" include("InterfaceI/ode_saveat_tests.jl")
    @time @safetestset "save_idxs Tests" include("InterfaceI/ode_saveidxs_tests.jl")
    @time @safetestset "Scalar Handling Tests" include("InterfaceI/scalar_handling_tests.jl")
    @time @safetestset "Static Array Tests" include("InterfaceI/static_array_tests.jl")
    @time @safetestset "derivative_discontinuity Tests" include("InterfaceI/derivative_discontinuity_test.jl")
    @time @safetestset "Composite Algorithm Tests" include("InterfaceI/composite_algorithm_test.jl")
    @time @safetestset "Complex Tests" include("InterfaceI/complex_tests.jl")
    @time @safetestset "Ndim Complex Tests" include("InterfaceI/ode_ndim_complex_tests.jl")
    @time @safetestset "Number Type Tests" include("InterfaceI/ode_numbertype_tests.jl")
    @time @safetestset "Interpolation Output Type Tests" include("InterfaceI/interpolation_output_types.jl")
    @time @safetestset "Stiffness Detection Tests" include("InterfaceI/stiffness_detection_test.jl")
    @time @safetestset "Composite Interpolation Tests" include("InterfaceI/composite_interpolation.jl")
    @time @safetestset "Export tests" include("InterfaceI/export_tests.jl")
    @time @safetestset "Type Handling Tests" include("InterfaceI/type_handling.jl")
    @time @safetestset "Controller Tests" include("InterfaceI/controllers.jl")
    @time @safetestset "Inplace Interpolation Tests" include("InterfaceI/inplace_interpolation.jl")
    @time @safetestset "Algebraic Interpolation Tests" include("InterfaceI/algebraic_interpolation.jl")
    @time @safetestset "Interpolation and Cache Stripping Tests" include("InterfaceI/ode_strip_test.jl")
    @time @safetestset "Public API Package Splits" include("InterfaceI/public_api_package_split.jl")
    @time @safetestset "Aliasing Tests" include("InterfaceI/aliasing_tests.jl")
    return @time @safetestset "Solution Memory Release" include("InterfaceI/solution_memory_tests.jl")
end

function interface_ii()
    is_APPVEYOR && return
    #@time @safetestset "No Recompile Tests" include("shared/norecompile.jl") # doesn't work on CI?
    @time @safetestset "AutoSparse Detection Tests" include("InterfaceII/autosparse_detection_tests.jl")
    @time @safetestset "Enum Tests" include("InterfaceII/enums.jl")
    return @time @safetestset "Get du Tests" include("InterfaceII/get_du.jl")
end

function interface_iii()
    is_APPVEYOR && return
    @time @safetestset "Derivative Utilities Tests" include("InterfaceIII/utility_tests.jl")
    @time @safetestset "stats Tests" include("InterfaceIII/stats_tests.jl")
    @time @safetestset "No Index Tests" include("InterfaceIII/noindex_tests.jl")
    @time @safetestset "Events + DAE addsteps Tests" include("InterfaceIII/event_dae_addsteps.jl")
    @time @safetestset "Units Tests" include("InterfaceIII/units_tests.jl")
    return @time @safetestset "DEVerbosity Tests" include("InterfaceIII/verbosity.jl")
end

function interface_iv()
    is_APPVEYOR && return
    @time @safetestset "Ambiguity Tests" include("InterfaceIV/ambiguity_tests.jl")
    @time @safetestset "Precision Mixing Tests" include("InterfaceIV/precision_mixing.jl")
    @time @safetestset "Sized Matrix Tests" include("InterfaceIV/sized_matrix_tests.jl")
    return @time @safetestset "Second Order with First Order Solver Tests" include("InterfaceIV/second_order_with_first_order_solvers.jl")
end

function interface_v()
    is_APPVEYOR && return
    @time @safetestset "Interpolation Derivative Error Tests" include("InterfaceV/interpolation_derivative_error_tests.jl")
    return @time @safetestset "GPU AutoDiff Interface Tests" include("InterfaceV/gpu_autodiff_interface_tests.jl")
end

function integrators_i()
    is_APPVEYOR && return
    @time @safetestset "Reinit Tests" include("Integrators_I/reinit_test.jl")
    @time @safetestset "Events Tests" include("Integrators_I/ode_event_tests.jl")
    @time @safetestset "Alg Events Tests" include("Integrators_I/alg_events_tests.jl")
    @time @safetestset "Discrete Callback Dual Tests" include("Integrators_I/discrete_callback_dual_test.jl")
    @time @safetestset "Callback Allocation Tests" include("Integrators_I/callback_allocation_tests.jl")
    @time @safetestset "Iterator Tests" include("Integrators_I/iterator_tests.jl")
    @time @safetestset "Integrator Interface Tests" include("Integrators_I/integrator_interface_tests.jl")
    @time @safetestset "Error Check Tests" include("Integrators_I/check_error.jl")
    @time @safetestset "Event Detection Tests" include("Integrators_I/event_detection_tests.jl")
    @time @safetestset "Event Repetition Detection Tests" include("Integrators_I/event_repeat_tests.jl")
    @time @safetestset "Multi-VCC Mask Tests" include("Integrators_I/multi_vcc_mask_tests.jl")
    return @time @safetestset "Step Limiter Tests" include("Integrators_I/step_limiter_test.jl")
end

function integrators_ii()
    is_APPVEYOR && return
    @time @safetestset "Integrator RNG Tests" include("Integrators_II/integrator_rng_tests.jl")
    @time @safetestset "Reverse Directioned Event Tests" include("Integrators_II/rev_events_tests.jl")
    @time @safetestset "Differentiation Direction Tests" include("Integrators_II/diffdir_tests.jl")
    @time @safetestset "Resize Tests" include("Integrators_II/resize_tests.jl")
    @time @safetestset "Cache Tests" include("Integrators_II/ode_cache_tests.jl")
    @time @safetestset "Add Steps Tests" include("Integrators_II/ode_add_steps_tests.jl")
    return @time @safetestset "IMEX Split Function Tests" include("Integrators_II/split_ode_tests.jl")
end

function regression_i()
    is_APPVEYOR && return
    @time @safetestset "Dense Tests" include("Regression_I/ode_dense_tests.jl")
    @time @safetestset "Special Interp Tests" include("Regression_I/special_interps.jl")
    @time @safetestset "Inplace Tests" include("Regression_I/ode_inplace_tests.jl")
    @time @safetestset "Adaptive Tests" include("Regression_I/ode_adaptive_tests.jl")
    return @time @safetestset "Hard DAE Tests" include("Regression_I/hard_dae.jl")
end

function regression_ii()
    is_APPVEYOR && return
    @time @safetestset "PSOS Energy Conservation Tests" include("Regression_II/psos_and_energy_conservation.jl")
    @time @safetestset "Unrolled Tests" include("Regression_II/ode_unrolled_comparison_tests.jl")
    @time @safetestset "IIP vs OOP Tests" include("Regression_II/iipvsoop_tests.jl")
    return @time @safetestset "Inference Tests" include("Regression_II/inference.jl")
end

function algconvergence_i()
    is_APPVEYOR && return
    return @time @safetestset "Non-autonomous Convergence Tests" include("AlgConvergence_I/non-autonomous_convergence_tests.jl")
end

function algconvergence_iii()
    is_APPVEYOR && return
    return @time @safetestset "Split Methods Tests" include("AlgConvergence_III/split_methods_tests.jl")
end

# AD / Downstream / ODEInterfaceRegression: use SciMLTesting.activate_group_env
# so Julia < 1.11 gets `develop_sources!` (the `[sources]` LTS backport already
# shipped in SciMLTesting). Developing only the monorepo root left registry
# Rosenbrock ≤2.3.2 paired with monorepo Differentiation 3.3.0 (unsatisfiable
# after the General retrocap).
function activate_downstream_env()
    return activate_group_env(
        joinpath(@__DIR__, "Downstream");
        parent = dirname(@__DIR__),
    )
end

function activate_odeinterface_env()
    return activate_group_env(
        joinpath(@__DIR__, "ODEInterfaceRegression");
        parent = dirname(@__DIR__),
    )
end

function activate_ad_env()
    return activate_group_env(
        joinpath(@__DIR__, "AD");
        parent = dirname(@__DIR__),
    )
end

function downstream_group()
    is_APPVEYOR && return
    activate_downstream_env()
    @time @safetestset "Measurements Tests" include("Downstream/measurements.jl")
    @time @safetestset "Time derivative Tests" include("Downstream/time_derivative_test.jl")
    return @time @safetestset "DynamicQuantities + Measurements Tests" include("Downstream/dynamicquantities_measurements.jl")
end

function ad_group()
    # Enzyme and Mooncake run on both AD lanes; Zygote runs on Julia <= 1.11.
    is_APPVEYOR && return
    activate_ad_env()
    @time @safetestset "AD Tests" include("AD/ad_tests.jl")
    @time @safetestset "Autodiff Events Tests" include("AD/autodiff_events.jl")
    return @time @safetestset "Discrete Adjoint Tests" include("AD/discrete_adjoints.jl")
end

function odeinterface_group()
    # Don't run ODEInterface tests on prerelease
    (is_APPVEYOR || !isempty(VERSION.prerelease)) && return
    activate_odeinterface_env()
    @time @safetestset "Init dt vs dorpri tests" include("ODEInterfaceRegression/init_dt_vs_dopri_tests.jl")
    return @time @safetestset "ODEInterface Regression Tests" include("ODEInterfaceRegression/odeinterface_regression.jl")
end

function qa_group()
    is_APPVEYOR && return
    @time @safetestset "Quality Assurance Tests" include("qa/qa_tests.jl")
    return @time @safetestset "Sublibrary Group Selection" include("qa/sublibrary_group_selection.jl")
end

function activate_gpu_env()
    Pkg.activate(joinpath(@__DIR__, "gpu"))
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

# Root (integration) GPU tests. The test/gpu environment pulls the full
# OrdinaryDiffEq stack plus the relevant sublibraries and runs on the self-hosted
# `gpu` runner (Julia 1, where [sources] is honored natively). Every file in
# test/gpu is included, so the set stays in sync as GPU tests are added or
# removed there without touching this dispatch.
function gpu_group()
    is_APPVEYOR && return
    activate_gpu_env()
    gpudir = joinpath(@__DIR__, "gpu")
    @testset "GPU Tests" begin
        for f in sort(filter(f -> endswith(f, ".jl"), readdir(gpudir)))
            @time @eval @safetestset $("GPU: " * f) include($(joinpath(gpudir, f)))
        end
    end
    return nothing
end

@time begin
    # Monorepo sublibrary routing. The root reads GROUP to pick a `lib/<sub>`
    # sublibrary, transitively develops its `[sources]` on Julia < 1.11, then
    # `Pkg.test`s it with the sub-group handed off via ODEDIFFEQ_TEST_GROUP.
    # This is kept as an explicit pre-step (rather than delegated to
    # `run_tests`'s built-in `lib_dir` path) so the sublibrary `Pkg.test`
    # invocation — `julia_args`, `force_latest_compatible_version = false`,
    # `allow_reresolve = true` — stays byte-for-byte identical to the previous
    # runtests.jl. (`run_tests`'s sublibrary path only passes `allow_reresolve`,
    # which would silently drop `--depwarn=yes`/`force_latest_compatible_version`.)
    # The root group dispatch (All / Interface / Integrators / Regression / QA /
    # functional groups) is delegated to `run_tests` below.
    base_group, test_group = detect_sublibrary_group(GROUP, LIB_DIR)

    if !isempty(base_group) && isdir(joinpath(LIB_DIR, base_group))
        Pkg.activate(joinpath(LIB_DIR, base_group))
        # On Julia < 1.11, the [sources] section in Project.toml is not supported.
        # Manually Pkg.develop local path dependencies so CI tests the PR branch code.
        # We resolve transitively: each developed dependency's own [sources] are also
        # developed, so that packages like OrdinaryDiffEqRosenbrockTableaus (a source
        # dependency of OrdinaryDiffEqRosenbrock) are correctly found even when testing
        # a higher-level sublibrary like OrdinaryDiffEqDefault.
        if VERSION < v"1.11.0-DEV.0"
            developed = Set{String}()
            # Never develop the active project: when sublibraries cyclically
            # reference each other via [sources] (e.g. DiffEqDevTools points
            # back at OrdinaryDiffEqCore), the transitive walk below would
            # otherwise try to `Pkg.develop` the active project itself, which
            # Pkg refuses with "package <X> has the same name or UUID as the
            # active project".
            push!(developed, normpath(joinpath(LIB_DIR, base_group)))
            specs = Pkg.PackageSpec[]
            queue = [joinpath(LIB_DIR, base_group)]
            while !isempty(queue)
                pkg_dir = popfirst!(queue)
                toml_path = joinpath(pkg_dir, "Project.toml")
                isfile(toml_path) || continue
                toml = Pkg.TOML.parsefile(toml_path)
                if haskey(toml, "sources")
                    for (dep_name, source_spec) in toml["sources"]
                        if source_spec isa Dict && haskey(source_spec, "path")
                            dep_path = normpath(joinpath(pkg_dir, source_spec["path"]))
                            if isdir(dep_path) && !(dep_path in developed)
                                push!(developed, dep_path)
                                @info "Queuing local source dependency" dep_name dep_path
                                push!(specs, Pkg.PackageSpec(path = dep_path))
                                # Queue this dependency so its own [sources] are also resolved.
                                push!(queue, dep_path)
                            end
                        end
                    end
                end
            end
            # Batch the develop call so Pkg resolves all path deps together;
            # calling it one-at-a-time would re-resolve the active project and
            # fail to find unregistered siblings.
            isempty(specs) || Pkg.develop(specs)
        end
        withenv("ODEDIFFEQ_TEST_GROUP" => test_group) do
            Pkg.test(base_group, julia_args = ["--check-bounds=auto", "--compiled-modules=yes", "--depwarn=yes"], force_latest_compatible_version = false, allow_reresolve = true)
        end
    else
        # Root-package group dispatch. `run_tests` owns the All / Interface /
        # Integrators / Regression / QA / functional-group routing that the old
        # hand-written `if GROUP == ...` ladder expressed. The group bodies
        # `@safetestset include` files from the per-group `test/<GroupKey>/`
        # folders (the v1.2 folder layout).
        run_tests(;
            # No root "Core" body: the previous runtests.jl had no Core branch.
            # A no-op core keeps GROUP=Core a no-op (matching the old behavior)
            # and is excluded from `all`.
            core = () -> nothing,
            groups = Dict(
                "InterfaceI" => interface_i,
                "InterfaceII" => interface_ii,
                "InterfaceIII" => interface_iii,
                "InterfaceIV" => interface_iv,
                "InterfaceV" => interface_v,
                "Integrators_I" => integrators_i,
                "Integrators_II" => integrators_ii,
                "Regression_I" => regression_i,
                "Regression_II" => regression_ii,
                "AlgConvergence_I" => algconvergence_i,
                "AlgConvergence_III" => algconvergence_iii,
                "Downstream" => downstream_group,
                "AD" => ad_group,
                "ODEInterfaceRegression" => odeinterface_group,
                "GPU" => gpu_group,
            ),
            # QA runs in the root test environment (no per-group Project.toml);
            # its body is the ExplicitImports testset, not the standard
            # Aqua/JET `run_qa`.
            qa = qa_group,
            # `All` runs exactly the Interface/Integrators/Regression functional
            # groups (in this order), matching the previous `GROUP == "All"`
            # branches. It deliberately EXCLUDES QA, AlgConvergence_*,
            # Downstream, AD, and ODEInterfaceRegression — none of which the old
            # `All` branch ran.
            all = [
                "InterfaceI", "InterfaceII", "InterfaceIII", "InterfaceIV", "InterfaceV",
                "Integrators_I", "Integrators_II",
                "Regression_I", "Regression_II",
            ],
            # Umbrella groups: one GROUP value triggers several functional
            # groups, reproducing the old `GROUP == "Interface"` /
            # `"Integrators"` / `"Regression"` branches.
            umbrellas = Dict(
                "Interface" => ["InterfaceI", "InterfaceII", "InterfaceIII", "InterfaceIV", "InterfaceV"],
                "Integrators" => ["Integrators_I", "Integrators_II"],
                "Regression" => ["Regression_I", "Regression_II"],
            ),
            # Monorepo sublibrary handoff var: the root reads GROUP, but each
            # sublibrary reads ODEDIFFEQ_TEST_GROUP for its sub-group. The
            # sublibrary Pkg.test is done explicitly above; sublib_env/lib_dir
            # are passed for completeness so the routing config is self-describing.
            sublib_env = "ODEDIFFEQ_TEST_GROUP",
            lib_dir = LIB_DIR,
        )
    end
end # @time
