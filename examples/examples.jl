using Revise
using StayOnTheRidge
using Symbolics

# general hypercube mapping
H(x; a = -1, b = 1) = a .+ (b .- a) .* x
H_inverse(x; a = -1, b = 1) = (x .- a) ./ (b .- a)
P(x; a = -1, b = 1) = min.(max.(x, a), b)

# example 1 - example from https://proceedings.mlr.press/v195/daskalakis23b.html (Appendix D)
function example1()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-2
    ϵ = 1e-2
    f = (x[1]-1/2)*(x[2]-1/2)
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end

# example 2 - normal form game expected utility function (2 players, 2 strategies)
function example2()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-2
    A = [4 -5; -5 100] # this cycle when γ = 1e-4, ϵ = 1e-3
    f = A[1,1]*x[1]*x[2]+A[1,2]*x[1]*(1-x[2])+A[2,1]*(1-x[1])*x[2]+A[2,2]*(1-x[1])*(1-x[2])
    s = Settings(f, x, n,min_coords, γ, ϵ)
    return s
end

# example 3 - function f2 from https://proceedings.mlr.press/v195/daskalakis23b.html (Appendix E)
function example3()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-4
    f = ifelse((x[1]^2 + x[2]^2)/2 <= 0, -x[1]*x[2] - 1/20*x[2]^2,
    ifelse((x[1]^2 + x[2]^2)/2 < 1, -x[1]*x[2] - 1/20*x[2]^2 + 2/20*(3*((x[1]^2 + x[2]^2)/2)^2 - 2*((x[1]^2 + x[2]^2)/2)^3)*x[2]^2,
        -x[1]*x[2] - 1/20*x[2]^2 + 2/20*x[2]^2))
    s = Settings(f, x, n, min_coords, γ, ϵ, H, H_inverse, P)
    return s
end

function example4() # https://link.springer.com/article/10.1007/s10589-019-00141-6 example 5.3 i)
    n = 4
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2]
    γ = 1e-10
    ϵ = 1e-10
    p = x[1]^2+x[2]^2-x[3]^2-x[4]^2
    q = x[1]+x[4]+1
    f = p/q
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end

function example5() # https://link.springer.com/article/10.1007/s10589-019-00141-6 example 5.3 ii)
    n = 6
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    p = sum(x[i]+x[n+1-i] for i in 1:3) + x[1]^2*x[5]^2-x[4]^2*x[2]^2 + x[1]^2*x[6]^2-x[4]^2*x[3]^2 + x[2]^2*x[6]^2 - x[5]^2*x[3]^2
    q = x[1]^2 + x[5]^2 + x[3]*x[6] + 1
    f = p/q
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end # X CYKLI SE

function example6() # https://arxiv.org/pdf/2109.04178.pdf example 1
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [2]
    γ = 1e-3
    ϵ = 1e-2
    f = 2*x[1]*x[2]^2-x[1]^2-x[2]
    s = Settings(f, x, n, min_coords, γ, ϵ, H, H_inverse, P)
    return s
end

function example7() # f1 from https://proceedings.mlr.press/v195/daskalakis23b.html (Appendix E)
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-4
    f = (4*x[1]^2-(x[2]-3*x[1]+x[1]^3/20)^2-x[2]^4/10) * exp(-(x[1]^2+x[2]^2)/100)
    s = Settings(f, x, n, min_coords, γ, ϵ, H, H_inverse, P)
    return s
end # DAVA JINY VYSLEDEK

### 
function example8() # example 6.3 i) from https://arxiv.org/pdf/1809.01218.pdf
    n = 6
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    f = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
    s = Settings(f, x, n, min_coords, γ, ϵ, H, H_inverse, P)
    return s
end

function example9() # example 6.3 ii) from https://arxiv.org/pdf/1809.01218.pdf
    n = 6
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    f = x[1]^2 + x[2]^2 + x[3]^2 - x[4]^2 - x[5]^2 - x[6]^2 + x[1]*x[5] - x[2]*x[4] + x[1]*x[6] - x[3]*x[4] + x[2]*x[6] - x[3]*x[5]
    s = Settings(f, x, n, min_coords, γ, ϵ, H, H_inverse, P)
    return s
end # DAVA JINY VYSLEDEK

function example10() # example 6.2 i) from https://arxiv.org/pdf/1809.01218.pdf
    n = 4
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2]
    γ = 1e-1
    ϵ = 1e-1
    f = (x[1] + x[2] + x[3] + x[4] + 1)^2 -4*(x[1]*x[2] + x[2]*x[3] + x[3]*x[4] + x[4] + x[1])
    s = Settings(f, x, n, min_coords, γ, ϵ)
    return s
end # ASSUMPTIOIN 1 VIOLATED ALE JE TO POLYNOM ???

# settings1 = example1();
# min_max1, trajectory1 = run_dynamics(settings1)
# plot_trajectory2D(min_max1, trajectory1)

# settings2 = example2();
# min_max2, trajectory2 = run_dynamics(settings2)
# plot_trajectory2D(min_max2, trajectory2)

# settings3 = example3();
# min_max3, trajectory3 = run_dynamics(settings3)
# plot_trajectory2D(H(min_max3), H.(trajectory3))

# settings4 = example4();
# min_max4, trajectory4 = run_dynamics(settings4)
# plot_trajectory2D(min_max4, trajectory4)

# settings5 = example5();
# min_max5, trajectory5 = run_dynamics(settings5)
# plot_trajectory2D(min_max5, trajectory5)

# settings6 = example6();
# min_max6, trajectory6 = run_dynamics(settings6)
# plot_trajectory2D(H(min_max6), H.(trajectory6))

settings7 = example7();
min_max7, trajectory7 = run_dynamics(settings7)
plot_trajectory2D(H(min_max7), H.(trajectory7))

# settings8 = example8();
# min_max8, trajectory8 = run_dynamics(settings8)
# println(H(min_max8))

# settings9 = example9();
# min_max9, trajectory9 = run_dynamics(settings9)
# println(H(min_max9))

# settings10 = example10();
# min_max10, trajectory10 = run_dynamics(settings10)
# println(min_max10)
