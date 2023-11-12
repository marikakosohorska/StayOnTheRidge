using Revise
using Stonr
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
    min_max, trajectory = execute_stonr(s)
    return min_max, trajectory
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
    min_max, trajectory = execute_stonr(s)
    return min_max, trajectory
end

# example 3 - f2 from https://proceedings.mlr.press/v195/daskalakis23b.html (Appendix E)
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
    min_max, trajectory = execute_stonr(s)
    return min_max, trajectory
end

function example4() # https://link.springer.com/article/10.1007/s10589-019-00141-6 example 5.3 i)
    n = 4
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2]
    γ = 1e-2
    ϵ = 1e-2
    p = x[1]^2+x[2]^2-x[3]^2-x[4]^2
    q = x[1]+x[4]+1
    f = p/q
    s = Settings(f,x,n,min_coords,γ,ϵ)
    min_max, trajectory = execute_stonr(s)
    return min_max, trajectory
end

function example5() # https://link.springer.com/article/10.1007/s10589-019-00141-6 example 5.3 ii)
    n = 6
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2,3]
    γ = 1e-3
    ϵ = 1e-1
    p = sum(x[i]+x[n+1-i] for i in 1:3) + x[1]^2*x[5]^2-x[4]^2*x[2]^2 + x[1]^2*x[6]^2-x[4]^2*x[3]^2 + x[2]^2*x[6]^2 - x[5]^2*x[3]^2
    q = x[1]^2 + x[5]^2 + x[3]*x[6] + 1
    f = p/q
    s = Settings(f,x,n,min_coords,γ,ϵ)
    min_max, trajectory = execute_stonr(s)
    return min_max, trajectory
end

function example6() # https://arxiv.org/pdf/2109.04178.pdf example 1
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [2]
    γ = 1e-3
    ϵ = 1e-3
    f = 2*x[1]*x[2]^2-x[1]^2-x[2]
    s = Settings(f, x, n, min_coords, γ, ϵ, H, H_inverse, P)
    min_max, trajectory = execute_stonr(s)
    return min_max, trajectory
end

function example7() # f1 from https://proceedings.mlr.press/v195/daskalakis23b.html (Appendix E)
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-3
    ϵ = 1e-5
    f = (4*x[1]^2-(x[2]-3*x[1]+x[1]^3/20)^2-x[2]^4/10)*exp(-(x[1]^2+x[2]^2)/100)
    s = Settings(f, x, n, min_coords, γ, ϵ, H, H_inverse, P)
    min_max, trajectory = execute_stonr(s)
    return min_max, trajectory
end

function example8() # example 6.3 from https://arxiv.org/pdf/1809.01218.pdf
    n = 6
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1,2,3]
    γ = 1e-2
    ϵ = 1e-1
    f = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
    println(f)
    s = Settings(f, x, n, min_coords, γ, ϵ, H, H_inverse, P)
    min_max, trajectory = execute_stonr(s)
    return min_max, trajectory
end

#min_max, trajectory = example1()

# min_max, trajectory = example2()

# min_max, trajectory = example3()
# min_max = H(min_max)
# trajectory = H.(trajectory)

# min_max, trajectory = example4()

# min_max, trajectory = example5()

min_max, trajectory = example6()
min_max = H(min_max)
trajectory = H.(trajectory)

# min_max, trajectory = example7()
# min_max = H(min_max)
# trajectory = H.(trajectory)

# min_max, trajectory = example8()
# min_max = H(min_max)
# trajectory = H.(trajectory)

plot_trajectory2D(min_max, trajectory)