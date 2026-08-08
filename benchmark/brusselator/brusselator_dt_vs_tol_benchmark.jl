# Brusselator DT vs TOL Fair Comparison - Clean Version for PR
# Fixed timestep exponential (DT sweep) vs adaptive implicit (TOL sweep) for same total error
# Measures: err vs ref 1e-12, time clean after warmup, nf, true_matvecs via CountingMatrix
# This file uses matrix operator to avoid size(Function) bug we fixed in tests

using OrdinaryDiffEq, OrdinaryDiffEqBDF, OrdinaryDiffEqExponentialRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK
using LinearAlgebra, SparseArrays, SciMLBase, Printf
try using DiffEqBase: MatrixOperator catch; using SciMLOperators: MatrixOperator end
using ADTypes: AutoFiniteDiff

const NS = [16, 32, 64]
const IMPLICIT_TOLS = [1e-4, 1e-5, 1e-6]
const EXP_DTS = [0.05, 0.025, 0.0125, 0.005]
const KRYLOV_M = 30

const EXP_ALGS_DT = [
    ("ETDRK2_m30", () -> ETDRK2(krylov=true, m=KRYLOV_M)),
    ("ETDRK4_m30", () -> ETDRK4(krylov=true, m=KRYLOV_M)),
    ("HochOst4_m30", () -> HochOst4(krylov=true, m=KRYLOV_M)),
]

const ADAPTIVE_EXP_ALGS_TOL = [
    ("Exprb43_m30_fd", () -> Exprb43(autodiff=AutoFiniteDiff(), m=KRYLOV_M)),
]

const IMPLICIT_ALGS_TOL = [
    ("FBDF", () -> FBDF()),
    ("KenCarp4", () -> KenCarp4()),
]

function estimate_builds(name::String)
    if occursin("ETDRK2",name); 2 elseif occursin("ETDRK4",name); 5 elseif occursin("HochOst4",name); 9 else; 2 end
end

brusselator_f(x,y,t) = (((x-0.3)^2 + (y-0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
init_bruss(N::Int) = begin xyd=range(0,stop=1,length=N); u=zeros(N,N,2); for I in CartesianIndices((N,N)); x=xyd[I[1]]; y=xyd[I[2]]; u[I,1]=22*(y*(1-y))^(3/2); u[I,2]=27*(x*(1-x))^(3/2); end; u end

function build_brusselator_matrices(N::Int, alpha_scaled::Float64)
    n=N*N; L=spzeros(n,n)
    for j in 1:N, i in 1:N
        idx=(j-1)*N+i; ip1=i==N ? 1 : i+1; im1=i==1 ? N : i-1; jp1=j==N ? 1 : j+1; jm1=j==1 ? N : j-1
        idx_ip1=(j-1)*N+ip1; idx_im1=(j-1)*N+im1; idx_jp1=(jp1-1)*N+i; idx_jm1=(jm1-1)*N+i
        L[idx, idx]=-4*alpha_scaled; L[idx, idx_ip1]=alpha_scaled; L[idx, idx_im1]=alpha_scaled
        L[idx, idx_jp1]=alpha_scaled; L[idx, idx_jm1]=alpha_scaled
    end
    C=spzeros(2*n,2*n); C[1:n,1:n]=L; C[n+1:end,n+1:end]=L; return C
end

function bruss_nonlin_vec!(du,u,p,t)
    Nloc=p[5]; A_val=p[1]; B_val=p[2]
    u3=reshape(u,Nloc,Nloc,2); du3=reshape(du,Nloc,Nloc,2)
    xyd=range(0,stop=1,length=Nloc)
    @inbounds for II in CartesianIndices((Nloc,Nloc))
        i=II[1]; j=II[2]; x=xyd[i]; y=xyd[j]
        du3[i,j,1]=B_val+u3[i,j,1]^2*u3[i,j,2]-(A_val+1)*u3[i,j,1]+brusselator_f(x,y,t)
        du3[i,j,2]=A_val*u3[i,j,1]-u3[i,j,1]^2*u3[i,j,2]
    end
end

make_p(N::Int;A=3.4,B=1.0,α=10.0)=(A,B,α,1.0/N,N)
make_split_prob(N::Int; tspan=(0.0,1.0)) = begin p=make_p(N); αs=p[3]/p[4]^2; L=build_brusselator_matrices(N,αs); u0=vec(init_bruss(N)); f1=MatrixOperator(L); f2=(du,u,par,t)->bruss_nonlin_vec!(du,u,par,t); (SplitODEProblem(f1,f2,u0,tspan,p), L) end
make_full_prob(N::Int; tspan=(0.0,1.0)) = begin p=make_p(N); u0=vec(init_bruss(N)); function full!(du,u,par,t) Nloc=par[5]; α=par[3]; dx=par[4]; αs=α/dx^2; u3=reshape(u,Nloc,Nloc,2); du3=reshape(du,Nloc,Nloc,2); @inbounds for k in 1:2, j in 1:Nloc, i in 1:Nloc; ip1=i==Nloc ? 1 : i+1; im1=i==1 ? Nloc : i-1; jp1=j==Nloc ? 1 : j+1; jm1=j==1 ? Nloc : j-1; du3[i,j,k]=αs*(u3[im1,j,k]+u3[ip1,j,k]+u3[i,jp1,k]+u3[i,jm1,k]-4*u3[i,j,k]); end; tmp=similar(du); bruss_nonlin_vec!(tmp,u,par,t); du.+=tmp; end; ODEProblem(full!, u0, tspan, p) end

struct CountingMatrix <: AbstractMatrix{Float64}
    A::SparseMatrixCSC{Float64,Int}; count::Base.RefValue{Int}
end
CountingMatrix(A)=CountingMatrix(A,Ref(0))
Base.size(C::CountingMatrix)=size(C.A); Base.size(C::CountingMatrix,d::Int)=size(C.A,d)
Base.getindex(C::CountingMatrix,i::Int,j::Int)=getindex(C.A,i,j)
LinearAlgebra.mul!(y,C::CountingMatrix,x)=(C.count[]+=1; mul!(y,C.A,x))
LinearAlgebra.mul!(y,C::CountingMatrix,x,a,b)=(C.count[]+=1; mul!(y,C.A,x,a,b))
Base.:*(C::CountingMatrix,x::AbstractVector)=(C.count[]+=1; C.A*x)

function warmup_all(prob_full, prob_split)
    for (name,factory) in IMPLICIT_ALGS_TOL; try solve(prob_full, factory(), abstol=1e-4, reltol=1e-4, maxiters=10, save_everystep=false) catch; end end
    for (name,factory) in EXP_ALGS_DT; try solve(prob_split, factory(), dt=0.05, maxiters=10, save_everystep=false) catch; end end
end

function run_for_N(N::Int)
    @printf("\n==================== N=%d (%d ODEs) ====================\n", N, 2*N^2)
    prob_full=make_full_prob(N); prob_split,L=make_split_prob(N)
    t0=time(); sol_ref=solve(prob_full, FBDF(), abstol=1e-12, reltol=1e-12, maxiters=Int(2e6), save_everystep=false); t_ref=time()-t0
    @printf("Ref steps=%d time=%.2fs\n", sol_ref.stats.naccept, t_ref)
    u_ref=sol_ref.u[end]
    warmup_all(prob_full, prob_split)
    results=[]
    for tol in IMPLICIT_TOLS
        for (name,factory) in IMPLICIT_ALGS_TOL
            alg=factory()
            try
                t_s=time(); sol=solve(prob_full, alg, abstol=tol, reltol=tol, maxiters=Int(2e6), save_everystep=false); t_e=time()-t_s
                err=norm(sol.u[end]-u_ref)/max(1.0,norm(u_ref))
                @printf("  %-20s tol=%.0e err=%.2e time=%6.3f nf=%d acc=%d\n", name, tol, err, t_e, sol.stats.nf+sol.stats.nf2, sol.stats.naccept)
                push!(results,(N=N,param=tol,param_type="tol",alg=name,err=err,time=t_e,nf=sol.stats.nf+sol.stats.nf2,naccept=sol.stats.naccept,nreject=sol.stats.nreject,builds=0,m=0,est_matvecs=0,true_matvecs=0,dt=0.0))
            catch e; @printf("  %-20s tol=%.0e EX %s\n", name,tol,e); end
        end
        if N <= 16
            for (name,factory) in ADAPTIVE_EXP_ALGS_TOL
                alg=factory()
                try
                    t_s=time(); sol=solve(prob_full, alg, abstol=tol, reltol=tol, maxiters=Int(2e6), save_everystep=false); t_e=time()-t_s
                    err=norm(sol.u[end]-u_ref)/max(1.0,norm(u_ref))
                    @printf("  %-20s tol=%.0e err=%.2e time=%6.3f nf=%d acc=%d\n", name,tol,err,t_e,sol.stats.nf+sol.stats.nf2,sol.stats.naccept)
                    push!(results,(N=N,param=tol,param_type="tol",alg=name,err=err,time=t_e,nf=sol.stats.nf+sol.stats.nf2,naccept=sol.stats.naccept,nreject=sol.stats.nreject,builds=2,m=30,est_matvecs=2*30*sol.stats.naccept,true_matvecs=0,dt=0.0))
                catch e; @printf("  %-20s tol=%.0e EX %s\n", name,tol,e); end
            end
        end
    end
    for dt in [0.05, 0.025, 0.0125, 0.005]
        for (name,factory) in EXP_ALGS_DT
            alg=factory()
            try
                t_s=time(); sol=solve(prob_split, alg, dt=dt, maxiters=Int(2e6), save_everystep=false); t_e=time()-t_s
                err=norm(sol.u[end]-u_ref)/max(1.0,norm(u_ref))
                L_count=CountingMatrix(L); f1_count=MatrixOperator(L_count); u0=vec(init_bruss(N)); f2=(du,u,par,t)->bruss_nonlin_vec!(du,u,par,t); prob_count=SplitODEProblem(f1_count,f2,u0,(0.0,1.0),make_p(N)); L_count.count[]=0; sol_count=solve(prob_count, alg, dt=dt, maxiters=Int(2e6), save_everystep=false); true_mv=L_count.count[]
                @printf("  %-20s dt=%.4f err=%.2e time=%6.3f true_mv=%d\n", name, dt, err, t_e, true_mv)
                push!(results,(N=N,param=dt,param_type="dt",alg=name,err=err,time=t_e,nf=sol.stats.nf+sol.stats.nf2,naccept=sol.stats.naccept,nreject=sol.stats.nreject,builds=estimate_builds(name),m=30,est_matvecs=estimate_builds(name)*30*sol.stats.naccept,true_matvecs=true_mv,dt=dt))
            catch e; @printf("  %-20s dt=%.4f EX %s\n", name,dt,e); end
        end
    end
    return results
end

function estimate_builds(name::String)
    if occursin("ETDRK2",name); 2 elseif occursin("ETDRK4",name); 5 elseif occursin("HochOst4",name); 9 else; 2 end
end

all_results=[]; for N in [16,32,64]; append!(all_results, run_for_N(N)); end
csv_path=joinpath(@__DIR__, "brusselator_dt_vs_tol_with_adaptive_exp_warmup.csv")
open(csv_path,"w") do io; println(io,"N,param,param_type,alg,err,time,nf,naccept,nreject,builds,m,est_matvecs,true_matvecs,dt"); for r in all_results; println(io,"$(r.N),$(r.param),$(r.param_type),$(r.alg),$(r.err),$(r.time),$(r.nf),$(r.naccept),$(r.nreject),$(r.builds),$(r.m),$(r.est_matvecs),$(r.true_matvecs),$(r.dt)"); end; end
println("\nCSV saved: $csv_path")
