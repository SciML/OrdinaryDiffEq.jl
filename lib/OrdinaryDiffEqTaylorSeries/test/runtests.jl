using SciMLTesting
using OrdinaryDiffEqTaylorSeries, ODEProblemLibrary, DiffEqDevTools
using SciMLBase
using Test
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    using ForwardDiff
    using OrdinaryDiffEqCore: PIDController
    @time @safetestset "SciMLBase reexport" begin
        using OrdinaryDiffEqTaylorSeries, Test
        exported = (
            :ODEProblem, :ODEFunction, :SplitODEProblem, :solve, :init, :step!,
            :remake, :ReturnCode, :CallbackSet, :ContinuousCallback, :terminate!,
            :u_modified!, :add_tstop!, :get_du, :EnsembleProblem,
        )
        @test all(Base.isexported.(Ref(OrdinaryDiffEqTaylorSeries), exported))
        internal = (
            :build_solution, :isinplace, :has_jac, :AbstractODEProblem,
            :StandardODEProblem, :UJacobianWrapper, :LinearProblem,
            :ConvexOptimizationProblem,
        )
        @test !any(Base.isexported.(Ref(OrdinaryDiffEqTaylorSeries), internal))
    end
    @testset "ExplicitTaylor2 Convergence Tests" begin
        # Test convergence
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        sim = test_convergence(dts, prob_ode_linear, ExplicitTaylor2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob_ode_2Dlinear, ExplicitTaylor2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
    end

    # The in-place step evaluates `f` into a Taylor buffer whose value slot used to be
    # `k1` -- which is also the first-order partial of the input. An RHS that zeroes
    # `du` before accumulating into it therefore wiped its own Taylor seed, leaving
    # `k2 == 0` and dropping the method to first order. The library problems above
    # assign `du` directly, so they never exercised it.
    @testset "ExplicitTaylor2 with a zero-then-accumulate RHS" begin
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        analytic = (u0, p, t) -> u0 * exp(-t)

        direct! = (du, u, p, t) -> (@. du = -u; nothing)
        accumulate! = (du, u, p, t) -> (fill!(du, 0); @. du = du - u; nothing)

        for rhs! in (direct!, accumulate!)
            prob = ODEProblem(ODEFunction(rhs!; analytic = analytic), [1.0], (0.0, 1.0))
            sim = test_convergence(dts, prob, ExplicitTaylor2())
            @test sim.𝒪est[:final] ≈ 2 atol = testTol
        end

        # Both spellings are the same ODE, so they must give the same answer.
        pd = ODEProblem(ODEFunction(direct!; analytic = analytic), [1.0], (0.0, 1.0))
        pa = ODEProblem(ODEFunction(accumulate!; analytic = analytic), [1.0], (0.0, 1.0))
        sd = solve(pd, ExplicitTaylor2(); dt = 2.0^-6, adaptive = false)
        sa = solve(pa, ExplicitTaylor2(); dt = 2.0^-6, adaptive = false)
        @test sd.u[end] ≈ sa.u[end] rtol = 1.0e-12
    end

    @testset "ExplicitTaylorN Convergence Tests" begin
        # Test convergence
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        for N in 3:4
            alg = ExplicitTaylor(order = Val(N))
            sim = test_convergence(dts, prob_ode_linear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
            sim = test_convergence(dts, prob_ode_2Dlinear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
        end
    end

    @testset "ExplicitTaylorAdaptiveOrder Tests" begin
        sol = solve(prob_ode_linear, ExplicitTaylorAdaptiveOrder(min_order = Val(6), max_order = Val(10)))
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
        # a step count and a retcode alone are satisfied by a solver that never
        # leaves `u0`, which is exactly what issue #3976 was
        @test sol.u[end] ≈ prob_ode_linear.f.analytic(
            prob_ode_linear.u0, prob_ode_linear.p, sol.t[end]
        ) atol = 1.0e-6
    end

    # Issue #3976. Every order's jet writes into one Taylor buffer sized at
    # `max_order`. Storing a lower-order `TaylorScalar` there went through
    # TaylorDiff's number-promotion constructor, which keeps the constant term and
    # drops every derivative, so the step collapsed to `u = uprev`, the error
    # estimate was exactly zero and the solver reported Success after a handful of
    # enormous steps.
    @testset "ExplicitTaylorAdaptiveOrder keeps every Taylor coefficient" begin
        TaylorDiff = OrdinaryDiffEqTaylorSeries.TaylorDiff
        function f_exp!(du, u, p, t)
            du .= u
            return nothing
        end
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(f_exp!, [1.0, 2.0], (0.0, 1.0))
        integ = init(
            prob, ExplicitTaylorAdaptiveOrder(min_order = Val(2), max_order = Val(6)),
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        cache = integ.cache
        # `jets[i]` has order `min_order + i - 1`, and for u' = u it owes the buffer
        # that many nonzero coefficients on top of the value
        for (i, jet) in enumerate(cache.jets)
            fill!(cache.utaylor, zero(eltype(cache.utaylor)))
            jet(cache.utaylor, cache.coeffs[i], integ.uprev, integ.t)
            @test count(!iszero, TaylorDiff.flatten(cache.utaylor[1])) == i + 2
        end
    end

    @testset "ExplicitTaylorAdaptiveOrder accuracy tracks tolerance" begin
        for prob in (prob_ode_linear, prob_ode_2Dlinear)
            for tol in (1.0e-8, 1.0e-10)
                sol = solve(prob, ExplicitTaylorAdaptiveOrder(), abstol = tol, reltol = tol)
                @test SciMLBase.successful_retcode(sol)
                @test sol.u[end] != prob.u0
                exact = prob.f.analytic(prob.u0, prob.p, sol.t[end])
                @test maximum(abs.(sol.u[end] .- exact)) < 1.0e-6 * maximum(abs.(exact))
            end
        end

        # out-of-place array problems additionally hit an UndefRefError while
        # saving: `initialize!` sized `integrator.k` but `perform_step!` never
        # filled it
        prob_arr = ODEProblem((u, p, t) -> [-u[2], u[1]], [1.0, 0.0], (0.0, 1.0))
        sol = solve(
            prob_arr, ExplicitTaylorAdaptiveOrder(),
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end] ≈ [cos(1.0), sin(1.0)] atol = 1.0e-8
        @test sol(0.5) ≈ [cos(0.5), sin(0.5)] atol = 1.0e-8

        # the order window has to leave room for the extra jet the error estimate
        # rides on, otherwise the order loop never runs and a stale estimate is
        # what the step is accepted on
        @test_throws ArgumentError solve(
            prob_ode_linear,
            ExplicitTaylorAdaptiveOrder(min_order = Val(4), max_order = Val(4))
        )
    end

    # `test_convergence` solves with `adaptive = false`, which pins the order at
    # `max_order - 1` for the whole run. The jet order then equals the buffer
    # order, so this covers the top of the window only.
    @testset "ExplicitTaylorAdaptiveOrder Convergence Tests" begin
        dts = 2.0 .^ (-6:-2)
        testTol = 0.2
        for N in 3:5
            alg = ExplicitTaylorAdaptiveOrder(min_order = Val(N - 2), max_order = Val(N))
            sim = test_convergence(dts, prob_ode_linear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
            sim = test_convergence(dts, prob_ode_2Dlinear, alg)
            @test sim.𝒪est[:final] ≈ N atol = testTol
        end
    end

    # The mismatch of #3976 only exists below the top of the window, so hold the
    # order there by hand and check that the method still converges at the order
    # its jet claims. With the coefficients dropped the step is a no-op and the
    # error stops depending on `dt` at all.
    @testset "ExplicitTaylorAdaptiveOrder Convergence Below max_order" begin
        alg = ExplicitTaylorAdaptiveOrder(min_order = Val(2), max_order = Val(6))
        dts = 2.0 .^ (-8:-4)
        for prob in (prob_ode_linear, prob_ode_2Dlinear)
            exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
            errors = map(dts) do dt
                integrator = init(prob, alg; dt, adaptive = false)
                integrator.cache.current_order[] = 3
                solve!(integrator)
                maximum(abs.(integrator.sol.u[end] .- exact))
            end
            slopes = diff(log2.(errors)) ./ diff(log2.(dts))
            @test sum(slopes) / length(slopes) ≈ 4 atol = 0.2
        end
    end

    # Dense output on the in-place cache reads `max_order` coefficients back out
    # of `integrator.k`, which the adaptive-order step has to fill and the
    # interpolant has to read at the buffer order rather than at `min_order`.
    @testset "ExplicitTaylorAdaptiveOrder Dense Output" begin
        for prob in (prob_ode_linear, prob_ode_2Dlinear)
            sol = solve(
                prob, ExplicitTaylorAdaptiveOrder(),
                abstol = 1.0e-10, reltol = 1.0e-10, dense = true
            )
            for t in (0.25, 0.5, 0.75)
                exact = prob.f.analytic(prob.u0, prob.p, t)
                @test maximum(abs.(sol(t) .- exact)) < 1.0e-8 * maximum(abs.(exact))
            end
        end
    end

    # Issue #3976: once the order dropped off `max_order` the step froze, so the
    # solver walked to the end of the span in a handful of steps and reported
    # Success on a 96% error (7 steps here, 60 once the coefficients survive).
    @testset "ExplicitTaylorAdaptiveOrder on Pleiades" begin
        prob = remake(prob_ode_pleiades, tspan = (0.0, 1.0))
        ref = solve(
            prob, ExplicitTaylor(order = Val(6)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        sol = solve(
            prob, ExplicitTaylorAdaptiveOrder(min_order = Val(2), max_order = Val(6)),
            abstol = 1.0e-8, reltol = 1.0e-8
        )
        @test SciMLBase.successful_retcode(sol)
        @test length(sol.t) > 30
        @test maximum(abs.(sol.u[end] .- ref.u[end])) < 1.0e-6 * maximum(abs.(ref.u[end]))
    end

    # `get_fsalfirstlast` used to alias the FSAL buffer to `cache.u`, which is
    # `integrator.u`. The auto-dt heuristic writes an extra `f` evaluation into
    # `fsallast`, so the aliasing made that call `f(u, u, p, t)` and corrupted the
    # state before the first step, for any RHS that reads `u` after writing `du`.
    @testset "auto-dt must not write through the state buffer" begin
        function f_alias!(du, u, p, t)
            du[1] = 7.0
            du[2] = 9.0
            s = u[1] + u[2]
            du[1] = -s
            du[2] = s
            return nothing
        end
        u0 = [1.0, 2.0]
        for alg in (ExplicitTaylor(order = Val(4)), ExplicitTaylorAdaptiveOrder())
            integ = init(
                ODEProblem(f_alias!, copy(u0), (0.0, 1.0)), alg,
                abstol = 1.0e-8, reltol = 1.0e-8
            )
            @test !OrdinaryDiffEqTaylorSeries.isfsal(alg)
            @test integ.u == u0
        end
    end

    # End-to-end: with the state corrupted, the auto-dt heuristic saw a NaN and
    # collapsed `dt0` to `dtmin`, so the controller burned hundreds of steps
    # climbing back out (905 before the fix, 587 after, at this tolerance).
    @testset "Auto-dt on Pleiades does not collapse to dtmin" begin
        sol = solve(
            prob_ode_pleiades, ExplicitTaylor(order = Val(6)),
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        @test SciMLBase.successful_retcode(sol)
        @test all(isfinite, sol.u[end])
        @test length(sol.t) < 700
    end

    @testset "AdaptiveOrder Dual tspan derivative" begin
        function f_dual_t!(du, u, p, t)
            du[1] = -0.5 * u[1]
            du[2] = -1.5 * u[2]
            return nothing
        end
        function u1_at_tend(t1)
            prob = ODEProblem{true, SciMLBase.FullSpecialize}(
                f_dual_t!, [1.0, 1.0], (zero(t1), t1)
            )
            sol = solve(
                prob, ExplicitTaylorAdaptiveOrder(),
                abstol = 1.0e-8, reltol = 1.0e-8, save_everystep = false
            )
            return sol.u[end][1]
        end
        d = ForwardDiff.derivative(u1_at_tend, 1.0)
        @test isfinite(d)
        @test d ≈ -0.5 * exp(-0.5) rtol = 1.0e-6
    end

    # Recorded from master e70688847: Float32 t / Float64 u, PIDController(0.7, -0.4),
    # rates -0.7/-1.3, tol 1e-9. The tType-typed snapshot (15702e7e6) rounds QT=Float64
    # controller state and diverges by step 3 (68 vs 72 saved times).
    const MASTER_F32_PID_T = Float32[
        0.0f0, 0.063194774f0, 0.22562528f0, 0.33316094f0, 0.49529004f0, 0.5823058f0,
        0.8021585f0, 1.000312f0, 1.1055535f0, 1.2704147f0, 1.3592904f0, 1.5845374f0,
        1.8223977f0, 1.940033f0, 2.112692f0, 2.2070622f0, 2.4466908f0, 2.7284784f0,
        2.852737f0, 3.072989f0, 3.31804f0, 3.4523935f0, 3.6809213f0, 3.8024008f0,
        4.047207f0, 4.4841695f0, 5.0226393f0, 5.5167813f0, 5.7618885f0, 6.0359616f0,
        6.1867046f0, 6.563721f0, 6.88412f0, 7.0580525f0, 7.4938827f0, 7.847814f0,
        8.041755f0, 8.527738f0, 8.957954f0, 9.189331f0, 9.654674f0, 10.31106f0,
        10.657045f0, 11.241399f0, 11.779732f0, 12.069865f0, 12.799378f0, 13.405901f0,
        13.73935f0, 14.416192f0, 15.408678f0, 15.925745f0, 16.604174f0, 17.272423f0,
        18.009304f0, 18.806004f0, 19.631573f0, 20.496643f0, 21.425524f0, 22.424778f0,
        23.497795f0, 24.659286f0, 25.930174f0, 27.333923f0, 28.902752f0, 30.685595f0,
        32.758316f0, 35.24815f0, 37.197197f0, 42.084785f0, 48.60787f0, 50.0f0,
    ]
    const MASTER_F32_PID_U = [
        [1.0, 1.0],
        [0.9567278159155104, 0.9211308248226885],
        [0.853902988388895, 0.7457887842874269],
        [0.7919851477879968, 0.6484896801129594],
        [0.7070152714023238, 0.5252520563858398],
        [0.6652355788961215, 0.4690726980681235],
        [0.5703466720971561, 0.35246428924383333],
        [0.49647686899307447, 0.2724212842867606],
        [0.4612166172931189, 0.23758743643863955],
        [0.41094729765955373, 0.1917545456773377],
        [0.3861600845692751, 0.17083151529253418],
        [0.3298305646265314, 0.1274668878315592],
        [0.2792414486602964, 0.09356329241447084],
        [0.2571686695013546, 0.08029540243840652],
        [0.22789175084782373, 0.06415200974679519],
        [0.21332389961820403, 0.05674538229236388],
        [0.18038107747603715, 0.04155663611289715],
        [0.1480900188944324, 0.02881030745650024],
        [0.1357533085090612, 0.024512840666924902],
        [0.11635698415790129, 0.01840954831748047],
        [0.09801565974604877, 0.01338725643703886],
        [0.0892177262070979, 0.011241851683331354],
        [0.07602865396570299, 0.008352471238271797],
        [0.06983076548603741, 0.007132302863967219],
        [0.05883344310366911, 0.005188197412405759],
        [0.043329628032158986, 0.002939780647026137],
        [0.029722597952144044, 0.0014598357521645375],
        [0.021031224895453456, 0.0007679269531398939],
        [0.0177153789218489, 0.0005583885992665131],
        [0.014622801331834988, 0.00039102060900323137],
        [0.01315842286683989, 0.0003214347630942863],
        [0.010106213494514836, 0.00019689547261154597],
        [0.008075793590408123, 0.00012982070546088876],
        [0.007150042068210807, 0.00010354867361325272],
        [0.005270035927585881, 5.8760078887939035e-5],
        [0.004113545403781687, 3.7090172127025184e-5],
        [0.0035913451172941892, 2.8824599221587918e-5],
        [0.0025557328699633997, 1.5324485415788434e-5],
        [0.0018911546728389715, 8.759778830083749e-6],
        [0.0016083744162123965, 6.484282235318752e-6],
        [0.00116123495766494, 3.5410933650414185e-6],
        [0.0007334571721001263, 1.508527306269208e-6],
        [0.000575696005089054, 9.62090713896496e-7],
        [0.00038242482162196504, 4.5009053984970917e-7],
        [0.00026235521436837475, 2.2354585529805283e-7],
        [0.0002141348045320148, 1.533070380153576e-7],
        [0.0001285022765277209, 5.938726208985803e-8],
        [8.404738481529716e-5, 2.6993611443166422e-8],
        [6.655083044382362e-5, 1.749851017113533e-8],
        [4.143711664682195e-5, 7.25885210476706e-9],
        [2.068557437068707e-5, 1.9976885631789582e-9],
        [1.4403761763623945e-5, 1.020000583726836e-9],
        [8.958375011298813e-6, 4.2225108634880974e-10],
        [5.611476495877283e-6, 1.7712845187164188e-10],
        [3.35012487581904e-6, 6.79607785805188e-11],
        [1.9180491095100996e-6, 2.412437234011767e-11],
        [1.0761713250785822e-6, 8.248092385931682e-12],
        [5.873468582719857e-7, 2.6788576355676956e-12],
        [3.0655552423283245e-7, 8.00791908799808e-13],
        [1.523104662300837e-7, 2.1845328185650314e-13],
        [7.186643525592279e-8, 5.414437350219847e-14],
        [3.18731732612646e-8, 1.1961861375226012e-14],
        [1.3093868226671346e-8, 2.2924003252020745e-15],
        [4.901389829559685e-9, 3.696585843648956e-16],
        [1.6345023815644297e-9, 4.811138024149123e-17],
        [4.692352103065251e-10, 4.74943127018766e-18],
        [1.0996969583000887e-10, 3.2620432062118605e-19],
        [1.924741154209957e-11, 1.54364394973131e-20],
        [4.919466803649046e-12, 1.3705412781767754e-21],
        [-1.4851868073949582e-12, -5.190728587062748e-20],
        [-2.416345451746056e-12, -1.7175587452472575e-17],
        [-9.118459772333861e-13, -2.770806004198636e-18],
    ]

    @testset "AdaptiveOrder Float32 t PIDController matches master" begin
        function f_f32!(du, u, p, t)
            du[1] = -0.7 * u[1]
            du[2] = -1.3 * u[2]
            return nothing
        end
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            f_f32!, [1.0, 1.0], (0.0f0, 50.0f0)
        )
        sol = solve(
            prob, ExplicitTaylorAdaptiveOrder(),
            controller = PIDController(0.7, -0.4),
            abstol = 1.0e-9, reltol = 1.0e-9
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.t == MASTER_F32_PID_T
        @test sol.u == MASTER_F32_PID_U
    end

    # Recorded from master e70688847: BigFloat state/time, default PIController.
    # isbits-only snapshots drop q11/errold (QT=BigFloat) and diverge (11 vs 17 times).
    const MASTER_BF_T = BigFloat[
        big"0.0",
        big"0.04757851636964158346526902255253472163536834601067494967998382098382059101513839",
        big"0.0856753757148559201266582616504333327699536121702479516927844681444677673745329",
        big"0.161869094405284593449436739846230555039124144489393955718385762465762120093323",
        big"0.2797441413640230038742020158377061422542685097624648123995954391474839623255956",
        big"0.4237953130912184030993616126701649138034148078099718499196763513179871698604765",
        big"0.4798052750715071184154935560131516330465421582797591153116557542091293645328977",
        big"0.5838442200572056947323948448105213064372337266238072922493246033673034916263216",
        big"0.6263969444979953948219895791287562779742205987367909195602932945638702534988932",
        big"0.7011175555983963158575916633030050223793764789894327356482346356510962027991163",
        big"0.730508972462021710184937269058022444242267801694738531778909560645582677273522",
        big"0.7844071779696173681479498152464179553938142413557908254231276485164661418188193",
        big"0.8228006864728077118150682276325465547626364266065770150168345331839046531176906",
        big"0.8582405491782461694520545821060372438101091412852616120134676265471237224941525",
        big"0.8865924393425969355616436656848297950480873130282092896107741012376989779953203",
        big"0.9092739514740775484493149325478638360384698504225674316886192809901591823962579",
        big"1.0",
    ]
    const MASTER_BF_U = Vector{BigFloat}[
        BigFloat[big"1.0", big"1.0"],
        BigFloat[big"0.9764914756614956873492630367472475814407227941000272575147956661872041953238339", big"0.9311193871313879380511334128963546360837056593060103606741327613474766023715165"],
        BigFloat[big"0.9580668833360280573870267200592879163381842680451089788975001399193943626751495", big"0.8794020742108145276902967264943756262200607904528028122693870250460325727933872"],
        BigFloat[big"0.9222540535768324029088944256053746018903798929839753276818690256190980070089457", big"0.7844255271853756580091200566092168549654676216873677639145441159400708192874692"],
        BigFloat[big"0.8694694589191709361208829985060902441659786191841463646298518183957248587010662", big"0.6572990348651466963766860644460164572451692022441926051400194289847517595527438"],
        BigFloat[big"0.8090474940448955437384521867091566318176815668592089991152208155728348868088575", big"0.529568386324605633684733600491631030464687355023854900723436416712210729625926"],
        BigFloat[big"0.7867044528220722802863008866299910914719807501421271612276392658451843454992704", big"0.4868944509229560523134423682800329116930514164327163598365493494274743588212096"],
        BigFloat[big"0.7468267040185886733131416593388451535771685977312108668132647232022724932259142", big"0.4165426881435485036199910453660386695344003490344004284052488605007038701102969"],
        BigFloat[big"0.7311047941556203827206232194991950274607263862882200096116392413065117404542839", big"0.3907859088136847836591197522098679228201575241112906962793768745534676585339942"],
        BigFloat[big"0.7042944356516384539269002394287998383187700620980001175170835805635439532975577", big"0.3493516281793388075386271739914164599627463197759497935887577089044866430306028"],
        BigFloat[big"0.694020009866022480807136740981932181804816531478569623309189932913609951785845", big"0.3342842972489301493729696285705194668945199073237676262510111460346409382866709"],
        BigFloat[big"0.6755665620485358776194124407677231102450496290408150297229881607914896459583509", big"0.3083219446916838282549422691596349293324663467484366475292674937377883470701665"],
        BigFloat[big"0.662721562389765535075425196793988614575728136204976401042624567943169728806761", big"0.2910672235547535053055116687354164713891661582814401514993203280711340003138237"],
        BigFloat[big"0.6510816158977401781777791262343281787404777689904087002648695069974431577846539", big"0.2759982307068389923168732117471268316141562111447460935554292915958037549341032"],
        BigFloat[big"0.6419170304703913993309454622749932546791046089370896219356683745187343644575798", big"0.264506710098130679189365230503424575394921183148869382648115194811218862936433"],
        BigFloat[big"0.6346783297651142843018009597098408092288625153396616987378662014168808876414782", big"0.2556589556388338666041176510956092599162079892962938741019327181968490359977441"],
        BigFloat[big"0.6065306597126333659558493742570283970116700448915072394836867434867447427363188", big"0.2231301601482997676260322732714902477638813900763077138378242466748202717823874"],
    ]

    @testset "AdaptiveOrder BigFloat PIController matches master" begin
        function f_bf!(du, u, p, t)
            du[1] = -0.5 * u[1]
            du[2] = -1.5 * u[2]
            return nothing
        end
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            f_bf!, [big(1.0), big(1.0)], (big(0.0), big(1.0))
        )
        sol = solve(
            prob, ExplicitTaylorAdaptiveOrder(),
            abstol = big(1.0e-10), reltol = big(1.0e-10)
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.t == MASTER_BF_T
        @test sol.u == MASTER_BF_U
    end

    @testset "AdaptiveOrder step! allocation bound" begin
        function f_alloc!(du, u, p, t)
            du[1] = -0.5 * u[1]
            du[2] = -1.5 * u[2]
            return nothing
        end
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            f_alloc!, [1.0, 1.0], (0.0, 1.0e6)
        )
        integrator = init(
            prob, ExplicitTaylorAdaptiveOrder(),
            abstol = 1.0e-8, reltol = 1.0e-8, save_everystep = false
        )
        for _ in 1:30
            step!(integrator)
        end
        @allocated step!(integrator)
        @allocated step!(integrator)
        allocs = @allocated step!(integrator)
        @test allocs < 2048
    end

    # Test AutoSpecialize (default ODEProblem wraps in FunctionWrappers)
    # and FullSpecialize paths for IIP problems
    @testset "AutoSpecialize / FullSpecialize IIP" begin
        # IIP array problem
        function f_iip!(du, u, p, t)
            du[1] = -u[2]
            du[2] = u[1]
            return nothing
        end
        u0 = [1.0, 0.0]
        tspan = (0.0, 1.0)

        # AutoSpecialize (default) - uses FunctionWrappers, unwrapped_f needed
        prob_auto = ODEProblem(f_iip!, u0, tspan)
        # FullSpecialize - no wrapping
        prob_full = ODEProblem{true, SciMLBase.FullSpecialize}(f_iip!, u0, tspan)

        for prob in (prob_auto, prob_full)
            sol2 = solve(prob, ExplicitTaylor2(), dt = 0.01)
            @test SciMLBase.successful_retcode(sol2)
            @test length(sol2.t) > 1

            sol8 = solve(
                prob, ExplicitTaylor(order = Val(8)),
                abstol = 1.0e-12, reltol = 1.0e-12
            )
            @test SciMLBase.successful_retcode(sol8)
            @test length(sol8.t) > 1
        end

        # Verify both give similar results
        sol_auto = solve(
            prob_auto, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        sol_full = solve(
            prob_full, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        @test sol_auto.u[end] ≈ sol_full.u[end] atol = 1.0e-10
    end

    # Test OOP (out-of-place) with array state
    @testset "OOP Array Problems" begin
        function f_oop(u, p, t)
            return [-u[2], u[1]]
        end
        u0 = [1.0, 0.0]
        tspan = (0.0, 1.0)

        prob_oop = ODEProblem(f_oop, u0, tspan)

        sol2 = solve(prob_oop, ExplicitTaylor2(), dt = 0.01)
        @test SciMLBase.successful_retcode(sol2)
        @test length(sol2.t) > 1

        sol8 = solve(
            prob_oop, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12
        )
        @test SciMLBase.successful_retcode(sol8)
        @test length(sol8.t) > 1

        # Check solution accuracy (harmonic oscillator: u1=cos(t), u2=sin(t))
        @test sol8.u[end][1] ≈ cos(1.0) atol = 1.0e-10
        @test sol8.u[end][2] ≈ sin(1.0) atol = 1.0e-10
    end

    # Dense output interpolation tests
    @testset "Taylor2 Dense Output" begin
        # Scalar OOP: u' = -u, u(0)=1 => u(t) = exp(-t)
        prob_scalar = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
        sol = solve(prob_scalar, ExplicitTaylor2(), dt = 0.01, dense = true)
        @test SciMLBase.successful_retcode(sol)
        # Check interpolation at intermediate points
        for t in 0.1:0.1:0.9
            @test sol(t) ≈ exp(-t) atol = 1.0e-3
        end

        # Array OOP: harmonic oscillator u1'=-u2, u2'=u1
        prob_arr = ODEProblem((u, p, t) -> [-u[2], u[1]], [1.0, 0.0], (0.0, 1.0))
        sol_arr = solve(prob_arr, ExplicitTaylor2(), dt = 0.01, dense = true)
        @test SciMLBase.successful_retcode(sol_arr)
        t_mid = 0.5
        @test sol_arr(t_mid)[1] ≈ cos(t_mid) atol = 1.0e-3
        @test sol_arr(t_mid)[2] ≈ sin(t_mid) atol = 1.0e-3

        # idxs parameter
        @test sol_arr(t_mid, idxs = 1) ≈ cos(t_mid) atol = 1.0e-3
        @test sol_arr(t_mid, idxs = 2) ≈ sin(t_mid) atol = 1.0e-3

        # IIP: same harmonic oscillator
        function f_dense_iip!(du, u, p, t)
            du[1] = -u[2]
            du[2] = u[1]
            return nothing
        end
        prob_iip = ODEProblem{true, SciMLBase.FullSpecialize}(
            f_dense_iip!, [1.0, 0.0], (0.0, 1.0)
        )
        sol_iip = solve(prob_iip, ExplicitTaylor2(), dt = 0.01, dense = true)
        @test SciMLBase.successful_retcode(sol_iip)
        @test sol_iip(t_mid)[1] ≈ cos(t_mid) atol = 1.0e-3
        @test sol_iip(t_mid)[2] ≈ sin(t_mid) atol = 1.0e-3
    end

    @testset "Taylor2 Dense Convergence" begin
        dts = 2.0 .^ (-8:-4)
        testTol = 0.2
        sim = test_convergence(
            dts, prob_ode_linear, ExplicitTaylor2(), dense_errors = true
        )
        @test sim.𝒪est[:L2] ≈ 2 atol = testTol
        sim = test_convergence(
            dts, prob_ode_2Dlinear, ExplicitTaylor2(), dense_errors = true
        )
        @test sim.𝒪est[:L2] ≈ 2 atol = testTol
    end

    @testset "TaylorN Dense Output" begin
        # Scalar OOP: u' = -u, u(0)=1 => u(t) = exp(-t)
        prob_scalar = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))

        for N in (4, 8)
            alg = ExplicitTaylor(order = Val(N))
            sol = solve(prob_scalar, alg, abstol = 1.0e-12, reltol = 1.0e-12, dense = true)
            @test SciMLBase.successful_retcode(sol)
            # High-order Taylor should give very accurate interpolation
            tol = N >= 8 ? 1.0e-10 : 1.0e-6
            for t in 0.1:0.1:0.9
                @test sol(t) ≈ exp(-t) atol = tol
            end
        end

        # Array OOP: harmonic oscillator
        prob_arr = ODEProblem((u, p, t) -> [-u[2], u[1]], [1.0, 0.0], (0.0, 1.0))
        sol8 = solve(
            prob_arr, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12, dense = true
        )
        @test SciMLBase.successful_retcode(sol8)
        t_mid = 0.5
        @test sol8(t_mid)[1] ≈ cos(t_mid) atol = 1.0e-10
        @test sol8(t_mid)[2] ≈ sin(t_mid) atol = 1.0e-10

        # idxs parameter
        @test sol8(t_mid, idxs = 1) ≈ cos(t_mid) atol = 1.0e-10
        @test sol8(t_mid, idxs = 2) ≈ sin(t_mid) atol = 1.0e-10

        # IIP: harmonic oscillator with AutoSpecialize
        function f_taylor_dense_iip!(du, u, p, t)
            du[1] = -u[2]
            du[2] = u[1]
            return nothing
        end
        prob_iip = ODEProblem(f_taylor_dense_iip!, [1.0, 0.0], (0.0, 1.0))
        sol_iip = solve(
            prob_iip, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-12, reltol = 1.0e-12, dense = true
        )
        @test SciMLBase.successful_retcode(sol_iip)
        @test sol_iip(t_mid)[1] ≈ cos(t_mid) atol = 1.0e-10
        @test sol_iip(t_mid)[2] ≈ sin(t_mid) atol = 1.0e-10

        # idxs on IIP
        @test sol_iip(t_mid, idxs = 1) ≈ cos(t_mid) atol = 1.0e-10
        @test sol_iip(t_mid, idxs = 2) ≈ sin(t_mid) atol = 1.0e-10
    end

    @testset "TaylorN Dense Convergence" begin
        dts = 2.0 .^ (-8:-4)
        testTol = 0.3
        for N in 3:4
            alg = ExplicitTaylor(order = Val(N))
            sim = test_convergence(
                dts, prob_ode_linear, alg, dense_errors = true
            )
            @test sim.𝒪est[:L2] ≈ N atol = testTol
            sim = test_convergence(
                dts, prob_ode_2Dlinear, alg, dense_errors = true
            )
            @test sim.𝒪est[:L2] ≈ N atol = testTol
        end
    end

    # Test IIP with a nonlinear system (Henon-Heiles-style)
    @testset "Nonlinear IIP System" begin
        function henon_heiles!(du, u, p, t)
            du[1] = -u[3] * (1 + 2u[4])
            du[2] = -u[4] - (u[3]^2 - u[4]^2)
            du[3] = u[1]
            du[4] = u[2]
            return nothing
        end
        u0 = [0.0, 0.5, 0.1, 0.0]
        tspan = (0.0, 10.0)

        prob = ODEProblem{true, SciMLBase.FullSpecialize}(henon_heiles!, u0, tspan)
        sol = solve(
            prob, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-14, reltol = 1.0e-14
        )
        @test SciMLBase.successful_retcode(sol)
        @test length(sol.t) > 1

        # Also test with AutoSpecialize
        prob_auto = ODEProblem(henon_heiles!, u0, tspan)
        sol_auto = solve(
            prob_auto, ExplicitTaylor(order = Val(8)),
            abstol = 1.0e-14, reltol = 1.0e-14
        )
        @test SciMLBase.successful_retcode(sol_auto)
        @test sol_auto.u[end] ≈ sol.u[end] atol = 1.0e-10
    end
end

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
