using DiffEqBase, Test, RecursiveArrayTools
import SciMLBase: AbstractSciMLFunction

macro iop_def(funcdef::Expr)
    """Define in- and out-of-place functions simultaneously.

    Call on oop function definition, defines two functions with suffixes _op and _ip.
    """
    @assert funcdef.head ∈ (:function, :(=)) && funcdef.args[1].head == :call

    fname = funcdef.args[1].args[1]

    opname = Symbol("$(fname)_op")
    ipname = Symbol("$(fname)_ip")

    opdef = deepcopy(funcdef)
    opdef.args[1].args[1] = opname

    return quote
        $(esc(opdef))
        $(esc(ipname))(du, args...) = du .= $(esc(opname))(args...)
    end
end

function test_inplace(du, expected, f::AbstractSciMLFunction, args...)
    """Test the in-place version of a function."""
    fill!(du, NaN)
    f(du, args...)
    return @test du == expected
end

# Allocate du automatically based on type of expected result
function test_inplace(expected, f::AbstractSciMLFunction, args...)
    return test_inplace(similar(expected), expected, f, args...)
end

function test_iop(expected, f_op::AbstractSciMLFunction, f_ip::AbstractSciMLFunction, args...)
    """Test in- and out-of-place version of function both match expected value."""
    @test f_op(args...) == expected
    return test_inplace(expected, f_ip, args...)
end

@iop_def f(u, p, t) = p[1] .* u
u = [1.0, 2.0, 3.0]
p = [2.0]
t = 0.0

@testset "ODEFunction with default recompile flag" begin
    odefun = ODEFunction{false}(f_op)
    odefun_ip = ODEFunction{true}(f_ip)
    expected = f_op(u, p, t)
    test_iop(expected, odefun, odefun_ip, u, p, t)
end

@testset "ODEFunction with recompile flag: $rflag" for rflag in (true, false)
    odefun = ODEFunction{false, rflag}(f_op)
    odefun_ip = ODEFunction{true, rflag}(f_ip)
    expected = f_op(u, p, t)
    test_iop(expected, odefun, odefun_ip, u, p, t)
end

# SplitFunction
@iop_def f2(u, p, t) = u .^ 2
sfun = SplitFunction{false}(f_op, f2_op)
sfun_ip = SplitFunction{true}(f_ip, f2_ip; _func_cache = similar(u))
expected = f_op(u, p, t) + f2_op(u, p, t)
test_iop(expected, sfun, sfun_ip, u, p, t)

# DynamicalODEFunction
@iop_def dode_f1(v, u, p, t) = -u
@iop_def dode_f2(v, u, p, t) = p[1] .* v
dodefun = DynamicalODEFunction{false}(dode_f1_op, dode_f2_op)
dodefun_ip = DynamicalODEFunction{true}(dode_f1_ip, dode_f2_ip)
v = [4.0, 5.0, 6.0]
expected = ArrayPartition(dode_f1_op(v, u, p, t), dode_f2_op(v, u, p, t))
test_iop(expected, dodefun, dodefun_ip, ArrayPartition(v, u), p, t)

# DiscreteFunction
dfun = DiscreteFunction{false}(f_op)
dfun_ip = DiscreteFunction{true}(f_ip)
test_iop(f_op(u, p, t), dfun, dfun_ip, u, p, t)

# Type stability
f_analytic(u, p, t) = u
jac = (u, p, t) -> 1
@inferred ODEFunction{false}(f_op, jac = jac)
@inferred DiscreteFunction{false}(f_op, analytic = f_analytic)
