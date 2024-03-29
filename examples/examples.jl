using Revise
using StayOnTheRidge
using Symbolics
using BenchmarkTools

# Example 1 - example from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix D)
# domain: [0,1]²
# solution given by this implementation: (0.51, 0.49)
# solution presented in the paper: (0.5, 0.5)
# characterization at [0,1]²: bilinear
function example1(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-2
    ϵ = 1e-2
    domain = Default()
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = (x[1]-1/2)*(x[2]-1/2)
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = (x[1]-1/2)*(x[2]-1/2)
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 2 - normal form game expected utility function (2 players, 2 strategies)
# domain: [0,1]²
# solution given by this implementation: (0.921, 0.921)
# solution computed by linear programming: (0.92105, 0.92105)
# characterization at [0,1]²: bilinear
function example2(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    A = [4 -5; -5 100]
    domain = Default()
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = A[1,1]*x[1]*x[2]+A[1,2]*x[1]*(1-x[2])+A[2,1]*(1-x[1])*x[2]+A[2,2]*(1-x[1])*(1-x[2])
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = A[1,1]*x[1]*x[2]+A[1,2]*x[1]*(1-x[2])+A[2,1]*(1-x[1])*x[2]+A[2,2]*(1-x[1])*(1-x[2])
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 3 - function f2 from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix E)
# domain: [-1,1]²
# solution given by this implementation: (-0.004, -0.002)
# solution presented in the paper: (0,0)
# characterization at [-1,1]²: neither convex in x₁ nor concave in x₂
function example3(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    domain = Hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = ifelse((x[1]^2 + x[2]^2)/2 <= 0, -x[1]*x[2] - 1/20*x[2]^2,
        ifelse((x[1]^2 + x[2]^2)/2 < 1, -x[1]*x[2] - 1/20*x[2]^2 + 2/20*(3*((x[1]^2 + x[2]^2)/2)^2 - 2*((x[1]^2 + x[2]^2)/2)^3)*x[2]^2,
            -x[1]*x[2] - 1/20*x[2]^2 + 2/20*x[2]^2))
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
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
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
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
    domain = Hyperrectangle([fill([0,1], 4)...])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        p = x[1]^2+x[2]^2-x[3]^2-x[4]^2
        q = x[1]+x[4]+1
        f = p/q
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        p_fd(x) = x[1]^2+x[2]^2-x[3]^2-x[4]^2
        q_fd(x) = x[1]+x[4]+1
        f_fd(x) = p_fd(x)/q_fd(x)
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 5 - https://arxiv.org/pdf/2109.04178.pdf example 1
# domain: [-1,1]²
# solution given by this implementation: (0.3967,0.6302)
# solution presented in the paper: (0.4,0.6)
# characterization at [-1,1]²: concave in x₁, neither convex nor concave in x₂
function example5(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    domain = Hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = -2*x[1]*x[2]^2+x[1]^2+x[2]
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = -2*x[1]*x[2]^2+x[1]^2+x[2]
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 6 - f1 from https://proceedings.mlr.press/v195/daskalakis23b/daskalakis23b.pdf (Appendix E)
# solution given by this implementation: (-1,-1)
# solution presented in the paper: (0,0)
# characterization at [-1,1]²: neither convex in x₁ nor concave in x₂
function example6(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-1
    domain = Hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x, i) for i in 1:n]
        f = (4*x[1]^2-(x[2]-3*x[1]+x[1]^3/20)^2-x[2]^4/10) * exp(-(x[1]^2+x[2]^2)/100)
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = (4*x[1]^2-(x[2]-3*x[1]+x[1]^3/20)^2-x[2]^4/10) * exp(-(x[1]^2+x[2]^2)/100)
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 7 - example 6.3 i) from https://arxiv.org/pdf/1809.01218.pdf
# domain: [-1,1]⁶
# solution given by this implementation: (-1,-1,1,1,1,1)
# solution presented in the paper: (-1,-1,1,1,1,1)
# characterization at [-1,1]⁶: neither convex in (x₁,x₂,x₃) nor concave in (x₄,x₅,x₆)
function example7(sym::Bool)
    n = 6
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    domain = Hyperrectangle([fill([-1,1], 6)...])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 8 - example 6.3 ii) from https://arxiv.org/pdf/1809.01218.pdf
# domain: [-1,1]⁶
# solution given by this implementation: (-0.0, -0.0014, -0.0013, -0.001, -0.0017, -0.0031)
# solution presented in the paper: (-1,1,-1,-1,1,-1)
# characterization at [-1,1]⁶: neither convex in (x₁,x₂,x₃) nor concave in (x₄,x₅,x₆)
function example8(sym::Bool)
    n = 6
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-2
    domain = Hyperrectangle([fill([-1,1], 6)...])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^2 + x[2]^2 + x[3]^2 - x[4]^2 - x[5]^2 - x[6]^2 + x[1]*x[5] - x[2]*x[4] + x[1]*x[6] - x[3]*x[4] + x[2]*x[6] - x[3]*x[5]
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = x[1]^2 + x[2]^2 + x[3]^2 - x[4]^2 - x[5]^2 - x[6]^2 + x[1]*x[5] - x[2]*x[4] + x[1]*x[6] - x[3]*x[4] + x[2]*x[6] - x[3]*x[5]
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 9 - function x₁² - x₂²
# domain: [-1,1]²
# solution given by this implementation: (-0.024, -0.024)
# well known solution: (0,0)
# characterization at [-1,1]²: convex in x₁, concave in x₂
function example9(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-1
    domain = Hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^2 - x[2]^2
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = x[1]^2 - x[2]^2
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 10 - monkey saddle
# domain: [-1,1]²
# solution given by this implementation: (-1,-1)
# well known solution: (0,0)
# characterization at [-1,1]²: neither convex in x₁ nor concave in x₂
function example10(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-1
    domain = Hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^3 - 3*x[1]*x[2]^2
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = x[1]^3 - 3*x[1]*x[2]^2
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 11 - example 1 from https://arxiv.org/pdf/2006.08141.pdf
# domain: [-1,1]x[-2π,2π]
# solution given by this implementation: (0,0)
# solutions presented in the paper: (0,-π), (0,π)
# characterization at [-1,1]x[-2π,2π]: convex in x₁, non-concave in x₂
function example11(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-4
    display(typeof([[-1,1],[-2*pi,2*pi]]))
    domain = Hyperrectangle([[-1,1],[-2*pi,2*pi]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = 0.2*x[1]*x[2]-cos(x[2])
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = 0.2*x[1]*x[2]-cos(x[2])
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 12 - example 3 from https://arxiv.org/pdf/2006.08141.pdf
# domain: [-1,1]²
# solution given by this implementation: (-1,1)
# solutions presented in the paper: (0,0), (-1,-1)
# characterization at [-1,1]²: non convex in x₁, concave in x₂
function example12(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-4
    domain = Hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = x[1]^3+2*x[1]*x[2]-x[2]^2
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = x[1]^3+2*x[1]*x[2]-x[2]^2
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 13 - figure 1 from https://arxiv.org/pdf/1807.02629.pdf
# domain: [0,1]²
# solution given by this implementation: (0.4189, 0.607)
# solutions presented in the paper: (0.4,0.6)
# characterization at [0,1]²: neither convex in x₁ nor concave in x₂
function example13(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    domain = Default()
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = (x[1]-0.5)*(x[2]-0.5)+1/3*exp(-(x[1]-1/4)^2-(x[2]-3/4)^2)
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = (x[1]-0.5)*(x[2]-0.5)+1/3*exp(-(x[1]-1/4)^2-(x[2]-3/4)^2)
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

# Example 14 - example 2.2 from https://arxiv.org/pdf/1807.02629.pdf
# domain: [-1,1]²
# solution given by this implementation: (-0.002,-0.002)
# solutions presented in the paper: (0,0)
# characterization at [-1,1]²: convex in x₁, non concave in x₂
function example14(sym::Bool)
    n = 2
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    domain = Hyperrectangle([[-1,1],[-1,1]])
    if sym
        x = [Symbolics.variable(:x,i) for i in 1:n]
        f = (x[1]^4*x[2]^2+x[1]^2+1)*(x[1]^2*x[2]^4-x[2]^2+1)
        return Config_sym(f, x, n, min_coords, γ, ϵ, domain)
    else
        f_fd(x) = (x[1]^4*x[2]^2+x[1]^2+1)*(x[1]^2*x[2]^4-x[2]^2+1)
        return Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
    end
end

config1 = example1(false);
elapsed1 = @elapsed min_max1, trajectory1, m1, k1 = run_dynamics(config1)
pretty_print(min_max1, elapsed1, m1, k1)
plot_trajectory2D(min_max1, trajectory1, config1.domain)
plot_contour2D(min_max1, config1.f, config1.domain)
plot_surface(min_max1, config1.f, config1.domain)

config2 = example2(false);
elapsed2 = @elapsed min_max2, trajectory2, m2, k2 = run_dynamics(config2)
pretty_print(min_max2, elapsed2, m2, k2)
plot_trajectory2D(min_max2, trajectory2, config2.domain)
plot_contour2D(min_max2, config2.f, config2.domain)
plot_surface(min_max2, config2.f, config2.domain)

config3 = example3(false);
elapsed3 = @elapsed min_max3, trajectory3, m3, k3 = run_dynamics(config3)
pretty_print(min_max3, elapsed3, m3, k3)
plot_trajectory2D(min_max3, trajectory3, config3.domain)
plot_contour2D(min_max3, config3.f, config3.domain)
plot_surface(min_max3, config3.f, config3.domain)

config4 = example4(false);
elapsed4 = @elapsed min_max4, trajectory4, m4, k4 = run_dynamics(config4)
pretty_print(min_max4, elapsed4, m4, k4)

config5 = example5(false);
elapsed5 = @elapsed min_max5, trajectory5, m5, k5 = run_dynamics(config5)
pretty_print(min_max5, elapsed5, m5, k5)
plot_trajectory2D(min_max5, trajectory5, config5.domain)
plot_contour2D(min_max5, config5.f, config5.domain)
plot_surface(min_max5, config5.f, config5.domain)

config6 = example6(false);
elapsed6 = @elapsed min_max6, trajectory6, m6, k6 = run_dynamics(config6)
pretty_print(min_max6, elapsed6, m6, k6)
plot_trajectory2D(min_max6, trajectory6, config6.domain)
plot_contour2D(min_max6, config6.f, config6.domain)
plot_surface(min_max6, config6.f, config6.domain)

config7 = example7(false);
elapsed7 = @elapsed min_max7, trajectory7, m7, k7 = run_dynamics(config7)
pretty_print(min_max7, elapsed7, m7, k7)

config8 = example8(false);
elapsed8 = @elapsed min_max8, trajectory8, m8, k8 = run_dynamics(config8)
pretty_print(min_max8, elapsed8, m8, k8)

config9 = example9(false);
elapsed9 = @elapsed min_max9, trajectory9, m9, k9 = run_dynamics(config9)
pretty_print(min_max9, elapsed9, m9, k9)
plot_trajectory2D(min_max9, trajectory9, config9.domain)
plot_contour2D(min_max9, config9.f, config9.domain)
plot_surface(min_max9, config9.f, config9.domain)

config10 = example10(false);
elapsed10 = @elapsed min_max10, trajectory10, m10, k10 = run_dynamics(config10)
pretty_print(min_max10, elapsed10, m10, k10)
plot_trajectory2D(min_max10, trajectory10, config10.domain)
plot_contour2D(min_max10, config10.f, config10.domain)
plot_surface(min_max10, config10.f, config10.domain)

config11 = example11(true);
elapsed11 = @elapsed min_max11, trajectory11, m11, k11 = run_dynamics(config11)
pretty_print(min_max11, elapsed11, m11, k11)
plot_trajectory2D(min_max11, trajectory11, config11.domain)
plot_contour2D(min_max11, config11.f, config11.domain)
plot_surface(min_max11, config11.f, config11.domain)

config12 = example12(false);
elapsed12 = @elapsed min_max12, trajectory12, m12, k12 = run_dynamics(config12)
pretty_print(min_max12, elapsed12, m12, k12)
plot_trajectory2D(min_max12, trajectory12, config12.domain)
plot_contour2D(min_max12, config12.f, config12.domain)
plot_surface(min_max12, config12.f, config12.domain)

config13 = example13(false);
elapsed13 = @elapsed min_max13, trajectory13, m13, k13 = run_dynamics(config13)
pretty_print(min_max13, elapsed13, m13, k13)
plot_trajectory2D(min_max13, trajectory13, config13.domain)
plot_contour2D(min_max13, config13.f, config13.domain)
plot_surface(min_max13, config13.f, config13.domain)

config14 = example14(false);
elapsed14 = @elapsed min_max14, trajectory14, m14, k14 = run_dynamics(config14)
pretty_print(min_max14, elapsed14, m14, k14)
plot_trajectory2D(min_max14, trajectory14, config14.domain)
plot_contour2D(min_max14, config14.f, config14.domain)
plot_surface(min_max14, config14.f, config14.domain)
