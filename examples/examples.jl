using Revise
using StayOnTheRidge
using Symbolics
using Plots
using LinearAlgebra
using BenchmarkTools

# general hypercube mapping
# function H_hypercube(a, b)
#     function H(x)
#         return a .+ (b .- a) .* x
#     end
# end

# general rectangle mapping
function H_hyperrectangle(sides)
    function H(x)
        return getindex.(sides,1) .+ (getindex.(sides,2) .- getindex.(sides, 1)) .* x
    end
end

# Example 1 - example from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix D)
# solution given by this implementation: (0.51, 0.49)
# solution presented in the paper: (0.5, 0.5)
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

# Example 2 - normal form game expected utility function (2 players, 2 strategies)
# solution given by this implementation: (0.921, 0.921)
# solution computed by linear programming: (0.92105, 0.92105)
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

# Example 3 - function f2 from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix E)
# solution given by this implementation: (-0.004, -0.002)
# solution presented in the paper: (0,0)
function example3(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    H = H_hyperrectangle([[-1,1],[-1,1]])
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

# Example 4 - https://link.springer.com/article/10.1007/s10589-019-00141-6 example 5.3 i)
# solution given by this implementation: (0,0,0,0)
# solution presented in the paper: (0.2540,0.2097,0.2487,0.2944)*10^(-4)
function example4(sym::Bool)
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

# Example 5 - https://arxiv.org/pdf/2109.04178.pdf example 1
# solution given by this implementation: (0.3967,0.6302)
# solution presented in the paper: (0.4,0.6)
function example5(sym::Bool)
    n = 2
    min_coords = [2]
    γ = 1e-4
    ϵ = 1e-3
    H = H_hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = 2*x[1]*x[2]^2-x[1]^2-x[2]
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = 2*x[1]*x[2]^2-x[1]^2-x[2]
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

# Example 6 - f1 from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix E)
# solution given by this implementation: (-1,-1)
# solution presented in the paper: (0,0)
function example6(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-1
    H = H_hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x, i) for i in 1:n]
        f = (4*x[1]^2-(x[2]-3*x[1]+x[1]^3/20)^2-x[2]^4/10) * exp(-(x[1]^2+x[2]^2)/100)
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = (4*x[1]^2-(x[2]-3*x[1]+x[1]^3/20)^2-x[2]^4/10) * exp(-(x[1]^2+x[2]^2)/100)
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

# Example 7 - example 6.3 i) from https://arxiv.org/pdf/1809.01218.pdf
# solution given by this implementation: (-1,-1,1,1,1,1)
# solution presented in the paper: (-1,-1,1,1,1,1)
function example7(sym::Bool)
    n = 6
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    H = H_hyperrectangle([fill([-1,1], 6)...])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

# Example 8 - example 6.3 ii) from https://arxiv.org/pdf/1809.01218.pdf
# solution given by this implementation: (0.0,-0.0141,-0.0142,-0.0107,-0.0177,-0.032)
# solution presented in the paper: (-1,1,-1,-1,1,-1)
function example8(sym::Bool)
    n = 6
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    H = H_hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^2 + x[2]^2 + x[3]^2 - x[4]^2 - x[5]^2 - x[6]^2 + x[1]*x[5] - x[2]*x[4] + x[1]*x[6] - x[3]*x[4] + x[2]*x[6] - x[3]*x[5]
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = x[1]^2 + x[2]^2 + x[3]^2 - x[4]^2 - x[5]^2 - x[6]^2 + x[1]*x[5] - x[2]*x[4] + x[1]*x[6] - x[3]*x[4] + x[2]*x[6] - x[3]*x[5]
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

# Example 9 - function x1^2 - x2^2
# solution given by this implementation: (-0.024, -0.024)
# well known solution: (0,0)
function example9(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-1
    H = H_hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^2 - x[2]^2
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = x[1]^2 - x[2]^2
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

# Example 10 - monkey saddle
# solution given by this implementation: (-1,-1)
# well known solution: (0,0)
function example10(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-1
    H = H_hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^3 - 3*x[1]*x[2]^2
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = x[1]^3 - 3*x[1]*x[2]^2
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

# Example 11 - example 1 from https://arxiv.org/pdf/2006.08141.pdf
# solution given by this implementation: (0,0)
# solutions presented in the paper: (0,-π), (0,π)
function example11(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-4
    H = H_hyperrectangle([[-1,1],[-2*pi,2*pi]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = 0.2*x[1]*x[2]-cos(x[2])
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = 0.2*x[1]*x[2]-cos(x[2])
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

# Example 12 - example 3 from https://arxiv.org/pdf/2006.08141.pdf
# solution given by this implementation: (-1,1)
# solutions presented in the paper: (0,0), (-1,-1)
function example12(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-4
    H = H_hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^3+2*x[1]*x[2]-x[2]^2
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = x[1]^3+2*x[1]*x[2]-x[2]^2
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

# Example 13 - figure 1 from https://arxiv.org/pdf/1807.02629.pdf
# solution given by this implementation: (0.4,0.6)
# solutions presented in the paper: (0.4,0.6)
function example13(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = (x[1]-0.5)*(x[2]-0.5)+1/3*exp(-(x[1]-1/4)^2-(x[2]-3/4)^2)
        return Config_sym(f, x, n, min_coords, γ, ϵ)
    else
        f_fd(x) = (x[1]-0.5)*(x[2]-0.5)+1/3*exp(-(x[1]-1/4)^2-(x[2]-3/4)^2)
        return Config_FD(f_fd, n, min_coords, γ, ϵ)
    end
end

# Example 14 - example 2.2 from https://arxiv.org/pdf/1807.02629.pdf
# solution given by this implementation: (-0.002,-0.002)
# solutions presented in the paper: (0,0)
function example14(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    H = H_hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = (x[1]^4*x[2]^2+x[1]^2+1)*(x[1]^2*x[2]^4-x[2]^2+1)
        return Config_sym(f, x, n, min_coords, γ, ϵ; H)
    else
        f_fd(x) = (x[1]^4*x[2]^2+x[1]^2+1)*(x[1]^2*x[2]^4-x[2]^2+1)
        return Config_FD(f_fd, n, min_coords, γ, ϵ; H)
    end
end

settings1 = example1(true);
elapsed1 = @elapsed min_max1, trajectory1, m1, k1 = run_dynamics(settings1)
pretty_print(min_max1, elapsed1, m1, k1)
plot_trajectory2D(min_max1, trajectory1, [0,1], [0,1])
plot_contour2D(min_max1, settings1.f, [0,1], [0,1])
plot_surface(min_max1, settings1.f, [0,1], [0,1])

settings2 = example2(true);
elapsed2 = @elapsed min_max2, trajectory2, m2, k2 = run_dynamics(settings2)
pretty_print(min_max2, elapsed2, m2, k2)
plot_trajectory2D(min_max2, trajectory2, [0,1], [0,1])
plot_contour2D(min_max2, settings2.f, [0,1], [0,1])
plot_surface(min_max2, settings2.f, [0,1], [0,1])

settings3 = example3(false);
elapsed3 = @elapsed min_max3, trajectory3, m3, k3 = run_dynamics(settings3)
pretty_print(settings3.H(min_max3), elapsed3, m3, k3)
plot_trajectory2D(settings3.H(min_max3), settings3.H.(trajectory3), [-1,1], [-1,1])
plot_contour2D(settings3.H(min_max3), settings3.f, [-1,1], [-1,1])
plot_surface(settings3.H(min_max3), settings3.f, [-1,1], [-1,1])

settings4 = example4(false);
elapsed4 = @elapsed min_max4, trajectory4, m4, k4 = run_dynamics(settings4)
pretty_print(min_max4, elapsed4, m4, k4)

settings5 = example5(false);
elapsed5 = @elapsed min_max5, trajectory5, m5, k5 = run_dynamics(settings5)
pretty_print(settings5.H(min_max5), elapsed5, m5, k5)
plot_trajectory2D(settings5.H(min_max5), settings5.H.(trajectory5), [-1,1], [-1,1])
plot_contour2D(settings5.H(min_max5), settings5.f, [-1,1], [-1,1])
plot_surface(settings5.H(min_max5), settings5.f, [-1,1], [-1,1])

settings6 = example6(true);
elapsed6 = @elapsed min_max6, trajectory6, m6, k6 = run_dynamics(settings6)
pretty_print(settings6.H(min_max6), elapsed6, m6, k6)
plot_trajectory2D(settings6.H(min_max6), settings6.H.(trajectory6), [-1,1], [-1,1])
plot_contour2D(settings6.H(min_max6), settings6.f, [-1,1], [-1,1])
plot_surface(settings6.H(min_max6), settings6.f, [-1,1], [-1,1])

settings7 = example7(false);
elapsed7 = @elapsed min_max7, trajectory7, m7, k7 = run_dynamics(settings7)
pretty_print(settings7.H(min_max7), elapsed7, m7, k7)

settings8 = example8(true);
elapsed8 = @elapsed min_max8, trajectory8, m8, k8 = run_dynamics(settings8)
pretty_print(settings8.H(min_max8), elapsed8, m8, k8)

settings9 = example9(true);
elapsed9 = @elapsed min_max9, trajectory9, m9, k9 = run_dynamics(settings9)
pretty_print(settings9.H(min_max9), elapsed9, m9, k9)
plot_trajectory2D(settings9.H(min_max9), settings9.H.(trajectory9), [-1,1], [-1,1])
plot_contour2D(settings9.H(min_max9), settings9.f, [-1,1], [-1,1])
plot_surface(settings9.H(min_max9), settings9.f, [-1,1], [-1,1])

settings10 = example10(true);
elapsed10 = @elapsed min_max10, trajectory10, m10, k10 = run_dynamics(settings10)
pretty_print(settings10.H(min_max10), elapsed10, m10, k10)
plot_trajectory2D(settings10.H(min_max10), settings10.H.(trajectory10), [-1,1], [-1,1])
plot_contour2D(settings10.H(min_max10), settings10.f, [-1,1], [-1,1])
plot_surface(settings10.H(min_max10), settings10.f, [-1,1], [-1,1])

settings11 = example11(true);
elapsed11 = @elapsed min_max11, trajectory11, m11, k11 = run_dynamics(settings11)
pretty_print(settings11.H(min_max11), elapsed11, m11, k11)
plot_trajectory2D(settings11.H(min_max11), settings11.H.(trajectory11), [-1,1], [-2*pi,2*pi])
plot_contour2D(settings11.H(min_max11), settings11.f, [-1,1], [-2*pi,2*pi])
plot_surface(settings11.H(min_max11), settings11.f, [-1,1], [-2*pi,2*pi])

settings12 = example12(true);
elapsed12 = @elapsed min_max12, trajectory12, m12, k12 = run_dynamics(settings12)
pretty_print(settings12.H(min_max12), elapsed12, m12, k12)
plot_trajectory2D(settings12.H(min_max12), settings12.H.(trajectory12), [-1,1], [-1,1])
plot_contour2D(settings12.H(min_max12), settings12.f, [-1,1], [-1,1])
plot_surface(settings12.H(min_max12), settings12.f, [-1,1], [-1,1])

settings13 = example13(true);
elapsed13 = @elapsed min_max13, trajectory13, m13, k13 = run_dynamics(settings13)
pretty_print(min_max13, elapsed13, m13, k13)
plot_trajectory2D(min_max13, trajectory13, [0,1], [0,1])
plot_contour2D(min_max13, settings13.f, [0,1], [0,1])
plot_surface(min_max13, settings13.f, [0,1], [0,1])

settings14 = example14(true);
elapsed14 = @elapsed min_max14, trajectory14, m14, k14 = run_dynamics(settings14)
pretty_print(settings14.H(min_max14), elapsed14, m14, k14)
plot_trajectory2D(settings14.H(min_max14), settings14.H.(trajectory14), [-1,1], [-1,1])
plot_contour2D(settings14.H(min_max14), settings14.f, [-1,1], [-1,1])
plot_surface(settings14.H(min_max14), settings14.f, [-1,1], [-1,1])


