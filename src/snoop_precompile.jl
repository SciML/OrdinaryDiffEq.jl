import SnoopPrecompile

SnoopPrecompile.@precompile_all_calls begin
    function lorenz(du, u, p, t)
        du[1] = 10.0(u[2] - u[1])
        du[2] = u[1] * (28.0 - u[3]) - u[2]
        du[3] = u[1] * u[2] - (8 / 3) * u[3]
    end

    function lorenz_oop(u, p, t)
        [10.0(u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2], u[1] * u[2] - (8 / 3) * u[3]]
    end

    solver_list = [
        BS3(), Tsit5(), Vern7(), Vern9(), Rosenbrock23(), Rosenbrock23(autodiff = false),
        #Rosenbrock23(chunk_size = 1), Rosenbrock23(chunk_size = Val{1}()),

        Rodas4(), Rodas4(autodiff = false),
        #Rodas4(chunk_size = 1), Rodas4(chunk_size = Val{1}()),

        Rodas5(), Rodas5(autodiff = false),
        #Rodas5(chunk_size = 1), Rodas5(chunk_size = Val{1}()),

        Rodas5P(), Rodas5P(autodiff = false),
        #Rodas5P(chunk_size = 1), Rodas5P(chunk_size = Val{1}()),

        TRBDF2(), TRBDF2(autodiff = false),
        #TRBDF2(chunk_size = 1), TRBDF2(chunk_size = Val{1}()),

        KenCarp4(), KenCarp4(autodiff = false),
        #KenCarp4(chunk_size = 1), KenCarp4(chunk_size = Val{1}()),

        QNDF(), QNDF(autodiff = false),
        #QNDF(chunk_size = 1), QNDF(chunk_size = Val{1}()),

        AutoTsit5(Rosenbrock23()), AutoTsit5(Rosenbrock23(autodiff = false)),
        #AutoTsit5(Rosenbrock23(chunk_size = 1)),
        #AutoTsit5(Rosenbrock23(chunk_size = Val{1}())),

        AutoTsit5(TRBDF2()), AutoTsit5(TRBDF2(autodiff = false)),
        #AutoTsit5(TRBDF2(chunk_size = 1)),
        #AutoTsit5(TRBDF2(chunk_size = Val{1}())),

        AutoVern9(KenCarp47()), AutoVern9(KenCarp47(autodiff = false)),
        #AutoVern9(KenCarp47(chunk_size = 1)),
        #AutoVern9(KenCarp47(chunk_size = Val{1}())),

        AutoVern9(Rodas5()), AutoVern9(Rodas5(autodiff = false)),
        #AutoVern9(Rodas5(chunk_size = 1)),
        #AutoVern9(Rodas5(chunk_size = Val{1}())),

        AutoVern9(Rodas5P()), AutoVern9(Rodas5P(autodiff = false)),
        #AutoVern9(Rodas5P(chunk_size = 1)),
        #AutoVern9(Rodas5P(chunk_size = Val{1}())),

        AutoVern7(Rodas4()), AutoVern7(Rodas4(autodiff = false)),
        #AutoVern7(Rodas4(chunk_size = 1)),
        #AutoVern7(Rodas4(chunk_size = Val{1}())),

        #AutoVern7(Rodas5P()), AutoVern7(Rodas5P(autodiff=false)),
        #AutoVern7(Rodas5P(chunk_size = 1)),
        #AutoVern7(Rodas5P(chunk_size = Val{1}())),

        AutoVern7(TRBDF2()), AutoVern7(TRBDF2(autodiff = false)),
        #AutoVern7(TRBDF2(chunk_size = 1)),
        #AutoVern7(TRBDF2(chunk_size = Val{1}())),
    ]

    prob_list = [ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                                                            (0.0, 1.0));
                 ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                                                            (0.0, 1.0), Float64[])
                 #ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0));
                 #ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[])
                 #ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0));
                 #ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[])
                 #ODEProblem{false,false}(lorenz_oop,[1.0;0.0;0.0],(0.0,1.0))
                 #ODEProblem{false,false}(lorenz_oop,[1.0;0.0;0.0],(0.0,1.0),Float64[])
                 ]

    for prob in prob_list, solver in solver_list
        solve(prob, solver)(5.0)
    end

    prob_list = nothing
end
