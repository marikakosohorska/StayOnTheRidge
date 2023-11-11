using Revise
using Stonr
using Symbolics

n = 2
x = [Symbolics.variable(:x,i) for i in 1:n]
min_coords = [1]
γ = 1e-4
ϵ = 1e-1

function determine_interval(point)
    return 1
end

function determine_interval_piecewise(point)
    val = (point[1]^2+point[2]^2)/2
    if val <= 0
        return 1
    elseif  val > 0 && val < 1
        return 2
    else
        return 3
    end
end

# example 1 - example from paper (Appendix D)
function example1()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-1
    cube_lb = 0
    cube_rb = 1
    f = (x[1]-1/2)*(x[2]-1/2)
    s = Settings([f],x,n,min_coords,γ,ϵ,cube_lb,cube_rb,determine_interval)
    min_max, trajectory = execute_stonr(s)
    plot_trajectory2D(min_max, trajectory,s)
end

# example 2 - general normal form game expected utility (2 players, 2 strategies)
function example2()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-1
    cube_lb = 0
    cube_rb = 1
    A = [4 -5; -5 100] # this cycle when γ = 1e-4, ϵ = 1e-3
    f = A[1,1]*x[1]*x[2]+A[1,2]*x[1]*(1-x[2])+A[2,1]*(1-x[1])*x[2]+A[2,2]*(1-x[1])*(1-x[2])
    println(f)
    s = Settings([f],x,n,min_coords,γ,ϵ,cube_lb,cube_rb,determine_interval)
    min_max, trajectory = execute_stonr(s)
    plot_trajectory2D(min_max, trajectory,s)
end

# example 3 - f2 from paper (Appendix E, piecewise function)
function example3()
    n = 2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    min_coords = [1]
    γ = 1e-4
    ϵ = 1e-1
    cube_lb = -1
    cube_rb = 1
    f1 = -x[1]*x[2]-1/20*x[2]^2
    f2 = -x[1]*x[2]-1/20*x[2]^2 + 2/20*(3*((x[1]^2+x[2]^2)/2)^2-2*((x[1]^2+x[2]^2)/2)^3)*x[2]^2
    f3 = -x[1]*x[2]-1/20*x[2]^2+2/20*x[2]^2
    s = Settings([f1,f2,f3],x,n,min_coords,γ,ϵ,cube_lb,cube_rb,determine_interval_piecewise)
    min_max, trajectory = execute_stonr(s)
    plot_trajectory2D(min_max, trajectory,s)
end

#example1()
#example2()
example3()