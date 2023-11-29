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

H_inverse(x; a = -1, b = 1) = (x .- a) ./ (b .- a)
P(x; a = -1, b = 1) = min.(max.(x, a), b)

# example 1 - example from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix D)
function example1(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-2
    ϵ = 1e-2
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = (x[1]-1/2)*(x[2]-1/2)
        return Config_sym(f, x, n, min_coords, γ, ϵ)
    else
        f_fd(x) = (x[1]-1/2)*(x[2]-1/2)
        return Config_FD(f_fd, n, min_coords, γ, ϵ)
    end
end

# example 2 - normal form game expected utility function (2 players, 2 strategies)
function example2(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    A = [4 -5; -5 100]
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = A[1,1]*x[1]*x[2]+A[1,2]*x[1]*(1-x[2])+A[2,1]*(1-x[1])*x[2]+A[2,2]*(1-x[1])*(1-x[2])
        return Config_sym(f, x, n, min_coords, γ, ϵ)
    else
        f_fd(x) = A[1,1]*x[1]*x[2]+A[1,2]*x[1]*(1-x[2])+A[2,1]*(1-x[1])*x[2]+A[2,2]*(1-x[1])*(1-x[2])
        return Config_FD(f_fd, n, min_coords, γ, ϵ)
    end
    return s
end

# example 3 - function f2 from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix E)
function example3(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    H = H_closure(-1,1)
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = ifelse((x[1]^2 + x[2]^2)/2 <= 0, -x[1]*x[2] - 1/20*x[2]^2,
        ifelse((x[1]^2 + x[2]^2)/2 < 1, -x[1]*x[2] - 1/20*x[2]^2 + 2/20*(3*((x[1]^2 + x[2]^2)/2)^2 - 2*((x[1]^2 + x[2]^2)/2)^3)*x[2]^2,
            -x[1]*x[2] - 1/20*x[2]^2 + 2/20*x[2]^2))
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        function f_fd(x)
            term1 = -x[1] * x[2] - 1/20 * x[2]^2
            term2 = -x[1] * x[2] - 1/20 * x[2]^2 + 2/20 * (3 * ((x[1]^2 + x[2]^2)/2)^2 - 2 * ((x[1]^2 + x[2]^2)/2)^3) * x[2]^2
            term3 = -x[1] * x[2] - 1/20 * x[2]^2 + 2/20 * x[2]^2
            if (x[1]^2 + x[2]^2)/2 <= 0
                return term1
            elseif (x[1]^2 + x[2]^2)/2 < 1
                return term2
            else
                return term3
            end
        end
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

function example4(sym::Bool) # https://link.springer.com/article/10.1007/s10589-019-00141-6 example 5.3 i)
    n = 4
    min_coords = [1,2]
    γ = 1e-3
    ϵ = 1e-2
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        p = x[1]^2+x[2]^2-x[3]^2-x[4]^2
        q = x[1]+x[4]+1
        f = p/q
        return Config_sym(f, x, n, min_coords, γ, ϵ)
    else
        p_fd(x) = x[1]^2+x[2]^2-x[3]^2-x[4]^2
        q_fd(x) = x[1]+x[4]+1
        f_fd(x) = p_fd(x)/q_fd(x)
        return Config_FD(f_fd, n, min_coords, γ, ϵ)
    end
end

function example5(sym::Bool) # https://arxiv.org/pdf/2109.04178.pdf example 1
    n = 2
    min_coords = [2]
    γ = 1e-4
    ϵ = 1e-3
    H = H_closure(-1,1)
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = 2*x[1]*x[2]^2-x[1]^2-x[2]
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = 2*x[1]*x[2]^2-x[1]^2-x[2]
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

function example6(sym::Bool) # f1 from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix E)
    n = 2
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-1
    H = H_closure(-1,1)
    if sym
        x = [Symbolics.variable(:x, i) for i in 1:n]
        f = (4*x[1]^2-(x[2]-3*x[1]+x[1]^3/20)^2-x[2]^4/10) * exp(-(x[1]^2+x[2]^2)/100)
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = (4*x[1]^2-(x[2]-3*x[1]+x[1]^3/20)^2-x[2]^4/10) * exp(-(x[1]^2+x[2]^2)/100)
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

function example7(sym::Bool) # example 6.3 i) from https://arxiv.org/pdf/1809.01218.pdf
    n = 6
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    H = H_closure(-1,1)
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

function example8(sym::Bool) # example 6.3 ii) from https://arxiv.org/pdf/1809.01218.pdf
    n = 6
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    H = H_closure(-1,1)
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^2 + x[2]^2 + x[3]^2 - x[4]^2 - x[5]^2 - x[6]^2 + x[1]*x[5] - x[2]*x[4] + x[1]*x[6] - x[3]*x[4] + x[2]*x[6] - x[3]*x[5]
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = x[1]^2 + x[2]^2 + x[3]^2 - x[4]^2 - x[5]^2 - x[6]^2 + x[1]*x[5] - x[2]*x[4] + x[1]*x[6] - x[3]*x[4] + x[2]*x[6] - x[3]*x[5]
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

function example9(sym::Bool) # function x1^2 - x2^2
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-1
    H = H_closure(-1,1)
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^2 - x[2]^2
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = x[1]^2 - x[2]^2
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

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

# settings1 = example1(false);
# elapsed1 = @elapsed min_max1, trajectory1, m1, k1 = run_dynamics(settings1)
# pretty_print(min_max1, elapsed1, m1, k1)
# plot_trajectory2D(min_max1, trajectory1, 0, 1)

# settings2 = example2(false);
# elapsed2 = @elapsed min_max2, trajectory2, m2, k2 = run_dynamics(settings2)
# pretty_print(min_max2, elapsed2, m2, k2)
# plot_trajectory2D(min_max2, trajectory2, 0, 1)

# settings3 = example3(false);
# elapsed3 = @elapsed min_max3, trajectory3, m3, k3 = run_dynamics(settings3)
# pretty_print(settings3.H(min_max3), elapsed3, m3, k3)
# plot_trajectory2D(settings3.H(min_max3), settings3.H.(trajectory3), -1, 1)

# settings4 = example4(false);
# elapsed4 = @elapsed min_max4, trajectory4, m4, k4 = run_dynamics(settings4)
# pretty_print(min_max4, elapsed4, m4, k4)
# plot_trajectory2D(min_max4, trajectory4, 0, 1)

settings5 = example5(false);
elapsed5 = @elapsed min_max5, trajectory5, m5, k5 = run_dynamics(settings5)
pretty_print(settings5.H(min_max5), elapsed5, m5, k5)
plot_trajectory2D(settings5.H(min_max5), settings5.H.(trajectory5), -1, 1)

# settings6 = example6(false);
# elapsed6 = @elapsed min_max6, trajectory6, m6, k6 = run_dynamics(settings6)
# pretty_print(settings6.H(min_max6), elapsed6, m6, k6)
# plot_trajectory2D(settings6.H(min_max6), settings6.H.(trajectory6), -1, 1)

# settings7 = example7(false);
# elapsed7 = @elapsed min_max7, trajectory7, m7, k7 = run_dynamics(settings7)
# pretty_print(settings7.H(min_max7), elapsed7, m7, k7)

# settings8 = example8(true);
# elapsed8 = @elapsed min_max8, trajectory8, m8, k8 = run_dynamics(settings8)
# pretty_print(settings8.H(min_max8), elapsed8, m8, k8)

# settings9 = example9(false);
# elapsed9 = @elapsed min_max9, trajectory9, m9, k9 = run_dynamics(settings9)
# pretty_print(settings9.H(min_max9), elapsed9, m9, k9)
# plot_trajectory2D(settings9.H(min_max9), settings9.H.(trajectory9), -1, 1)

