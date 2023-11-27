using Revise
using StayOnTheRidge
using Symbolics
using Plots
using LinearAlgebra
using BenchmarkTools

# general hypercube mapping
function H_closure(a, b)
    function H(x)
        return a .+ (b .- a) .* x
    end
end

#H(x; a = -1, b = 1) = a .+ (b .- a) .* x
H_inverse(x; a = -1, b = 1) = (x .- a) ./ (b .- a)
P(x; a = -1, b = 1) = min.(max.(x, a), b)

# example 1 - example from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix D)
function example1()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-2
    ϵ = 1e-2
    f = (x[1]-1/2)*(x[2]-1/2)
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end # ok

# example 2 - normal form game expected utility function (2 players, 2 strategies)
function example2()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    A = [4 -5; -5 100]
    f = A[1,1]*x[1]*x[2]+A[1,2]*x[1]*(1-x[2])+A[2,1]*(1-x[1])*x[2]+A[2,2]*(1-x[1])*(1-x[2])
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end # ok

# example 3 - function f2 from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix E)
function example3()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    f = ifelse((x[1]^2 + x[2]^2)/2 <= 0, -x[1]*x[2] - 1/20*x[2]^2,
    ifelse((x[1]^2 + x[2]^2)/2 < 1, -x[1]*x[2] - 1/20*x[2]^2 + 2/20*(3*((x[1]^2 + x[2]^2)/2)^2 - 2*((x[1]^2 + x[2]^2)/2)^3)*x[2]^2,
        -x[1]*x[2] - 1/20*x[2]^2 + 2/20*x[2]^2))
    H = H_closure(-1,1)
    s = Settings(f, x, n, min_coords, γ, ϵ; H)
    return s
end # ok

function example4() # https://link.springer.com/article/10.1007/s10589-019-00141-6 example 5.3 i)
    n = 4
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2]
    γ = 1e-3
    ϵ = 1e-2
    p = x[1]^2+x[2]^2-x[3]^2-x[4]^2
    q = x[1]+x[4]+1
    f = p/q
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end # ok

function example5() # https://link.springer.com/article/10.1007/s10589-019-00141-6 example 5.3 ii)
    n = 6
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2,3]
    γ = 1e-6
    ϵ = 1e-1
    p = sum(x[i]+x[n+1-i] for i in 1:3) + x[1]^2*x[5]^2-x[4]^2*x[2]^2 + x[1]^2*x[6]^2-x[4]^2*x[3]^2 + x[2]^2*x[6]^2 - x[5]^2*x[3]^2
    q = x[1]^2 + x[5]^2 + x[3]*x[6] + 1
    f = p/q
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end # cycle, no saddle point

function example6() # https://arxiv.org/pdf/2109.04178.pdf example 1
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [2]
    γ = 1e-4
    ϵ = 1e-3
    f = 2*x[1]*x[2]^2-x[1]^2-x[2]
    H = H_closure(-1,1)
    s = Settings(f, x, n, min_coords, γ, ϵ; H)
    return s
end # ok

function example7() # f1 from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix E)
    n = 2
    x = [Symbolics.variable(:x, i) for i in 1:n]
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-1
    f = (4*x[1]^2-(x[2]-3*x[1]+x[1]^3/20)^2-x[2]^4/10) * exp(-(x[1]^2+x[2]^2)/100)
    H = H_closure(-1,1)
    s = Settings(f, x, n, min_coords, γ, ϵ; H)
    return s
end # jiny vysledek

function example8() # example 6.3 i) from https://arxiv.org/pdf/1809.01218.pdf
    n = 6
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    f = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
    H = H_closure(-1,1)
    s = Settings(f, x, n, min_coords, γ, ϵ; H)
    return s
end # ok

function example9() # example 6.3 ii) from https://arxiv.org/pdf/1809.01218.pdf
    n = 6
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    f = x[1]^2 + x[2]^2 + x[3]^2 - x[4]^2 - x[5]^2 - x[6]^2 + x[1]*x[5] - x[2]*x[4] + x[1]*x[6] - x[3]*x[4] + x[2]*x[6] - x[3]*x[5]
    H = H_closure(-1,1)
    s = Settings(f, x, n, min_coords, γ, ϵ; H)
    return s
end # jiny vysledek, ale v reseni gradient 0 a hessian indefinitni

function example10() # example 6.2 i) from https://arxiv.org/pdf/1809.01218.pdf
    n = 4
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2]
    γ = 1e-5
    ϵ = 1e-1
    f = (x[1] + x[2] + x[3] + x[4] + 1)^2 -4*(x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4] + x[1])
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end # assumption 1 (and also 2) violated

function example11() # function x1^2 - x2^2
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-1
    f = x[1]^2 - x[2]^2
    H = H_closure(-1,1)
    s = Settings(f, x, n, min_coords, γ, ϵ; H)
    return s
end # ok

function example12()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-1
    f = -x[1]^2/2 + (x[2]+1)^2/2
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end # cycle

function example13() # example 5.1 from https://arxiv.org/pdf/2108.04698.pdf
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-2
    M = [0.8 -0.2; -0.2 0.5]
    R = sum(atan(x[i] - 5) for i in 1:2)
    f = 0.5*x'*M*x -5*R
    H = H_closure(-1,7)
    s = Settings(f, x, n, min_coords, γ, ϵ; H)
    display(f)
    return s
end # jiny vysledek

function format_elapsed(elapsed_seconds)
    if elapsed_seconds < 1e-6
        return string(round(elapsed_seconds * 1e9, digits=3), " ns")
    elseif elapsed_seconds < 1e-3
        return string(round(elapsed_seconds * 1e6, digits=3), " μs")
    elseif elapsed_seconds < 1
        return string(round(elapsed_seconds * 1e3, digits=3), " ms")
    else
        return string(round(elapsed_seconds, digits=3), " s")
    end
end

function pretty_print(point, time, m, k)
    str_point = "($(join(round.(point, digits=4), ", ")))"
    str_time = format_elapsed(time)
    display("Min max critical point: $str_point, elapsed time: $str_time, outer cycles: $m, inner cycles: $k")
end

settings1 = example1();
elapsed1 = @elapsed min_max1, trajectory1, m1, k1 = run_dynamics(settings1)
pretty_print(min_max1, elapsed1, m1, k1)
plot_trajectory2D(min_max1, trajectory1, 0, 1)

# settings2 = example2();
# elapsed2 = @elapsed min_max2, trajectory2, m2, k2 = run_dynamics(settings2)
# pretty_print(min_max2, elapsed2, m2, k2)
# plot_trajectory2D(min_max2, trajectory2, 0, 1)

# settings3 = example3();
# elapsed3 = @elapsed min_max3, trajectory3, m3, k3 = run_dynamics(settings3)
# pretty_print(settings3.H(min_max3), elapsed3, m3, k3)
# plot_trajectory2D(settings3.H(min_max3), settings3.H.(trajectory3), -1, 1)

# settings4 = example4();
# elapsed4 = @elapsed min_max4, trajectory4, m4, k4 = run_dynamics(settings4)
# pretty_print(min_max4, elapsed4, m4, k4)
# plot_trajectory2D(min_max4, trajectory4, 0, 1)

# settings5 = example5();
# elapsed5 = @elapsed min_max5, trajectory5, m5, k5 = run_dynamics(settings5)

# settings6 = example6();
# elapsed6 = @elapsed min_max6, trajectory6, m6, k6 = run_dynamics(settings6)
# pretty_print(settings6.H(min_max6), elapsed6, m6, k6)
# plot_trajectory2D(settings6.H(min_max6), settings6.H.(trajectory6), -1, 1)

# settings7 = example7();
# elapsed7 = @elapsed min_max7, trajectory7, m7, k7 = run_dynamics(settings7)
# pretty_print(settings7.H(min_max7), elapsed7, m7, k7)
# plot_trajectory2D(settings7.H(min_max7), settings7.H.(trajectory7), -1, 1)

# settings8 = example8();
# elapsed8 = @elapsed min_max8, trajectory8, m8, k8 = run_dynamics(settings8)
# pretty_print(settings8.H(min_max8), elapsed8, m8, k8)

# settings9 = example9();
# elapsed9 = @elapsed min_max9, trajectory9, m9, k9 = run_dynamics(settings9)
# pretty_print(settings9.H(min_max9), elapsed9, m9, k9)
# g = [settings9.grad_f[i](settings9.H(min_max9)) for i in 1:6]
# h = [settings9.hessian_f[i,j](settings9.H(min_max9)) for i in 1:6, j in 1:6]
# eigvals(h)

# settings10 = example10();
# elapsed10 = @elapsed min_max10, trajectory10, m10, k10 = run_dynamics(settings10)

# settings11 = example11();
# elapsed11 = @elapsed min_max11, trajectory11, m11, k11 = run_dynamics(settings11)
# pretty_print(settings11.H(min_max11), elapsed11, m11, k11)
# plot_trajectory2D(settings11.H(min_max11), settings11.H.(trajectory11), -1, 1)

# settings12 = example12();
# elapsed12 = @elapsed min_max12, trajectory12, m12, k12 = run_dynamics(settings12)

# settings13 = example13();
# elapsed13 = @elapsed min_max13, trajectory13, m13, k13 = run_dynamics(settings13)
# pretty_print(settings13.H(min_max13), elapsed13, m13, k13)
# plot_trajectory2D(settings13.H(min_max13), settings13.H.(trajectory13), -1, 7)

