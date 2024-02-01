using LinearAlgebra
using ForwardDiff
using Symbolics

include("domain.jl")
include("config.jl")

P(x; x_min = 0, x_max = 1) = min.(max.(x, x_min), x_max)

function get_V(conf::Config_sym, point)
    V = [conf.V[i](point) for i in 1:conf.n]
    return V
end

function get_V(conf::Config_FD, point)
    V = ForwardDiff.gradient(conf.f_obj, point)
    V[conf.min_coords] .*= -1
    return V
end

function get_Vi(conf::Config_sym, point, i)
    Vi = conf.V[i](point)
    return Vi
end

function get_Vi(conf::Config_FD, point, i)
    Vi = get_V(conf, point)[i]
    return Vi
end

function get_V_tr(conf::Config_sym, new_point, tr_coords)
    V_tr = [conf.V[tr](new_point) for tr in tr_coords]
    return V_tr
end

function get_V_tr(conf::Config_FD, new_point, tr_coords)
    V_tr = get_V(conf, new_point)[tr_coords]
    return V_tr
end

function get_jacobian_Si(conf::Config_sym, point, d_nonzero_idxs, S)
    jacobian_Si = [conf.J[i,j](point) for i in S, j in d_nonzero_idxs]
    return jacobian_Si
end

function get_jacobian_Si(conf::Config_FD, point, d_nonzero_idxs, S)
    jacobian = ForwardDiff.hessian(conf.f_obj, point)
    jacobian[conf.min_coords,:] .*= -1
    jacobian_Si = [jacobian[i,j] for i in S, j in d_nonzero_idxs]
    return jacobian_Si
end

function compute_direction(point, i, S, conf::Config)
    d = zeros(conf.n)

    if isempty(S)
        d[i] = 1
        return d
    end

    d_nonzero_idxs = vcat(S, [i])
    jacobian_Si = get_jacobian_Si(conf, point, d_nonzero_idxs, S)
    d_nonzero = nullspace(jacobian_Si)
    if size(d_nonzero, 2) != 1
        # direction is not uniquely defined
        error("Assumption 1 violated at i = $i, S = $S, x = $point, nullspace dimension is $(ndims(d_nonzero)), jacobian_Si_at_point = $jacobian_Si")
    end

    d_nonzero = vec(d_nonzero)
    decision_mat = hcat(transpose(jacobian_Si), d_nonzero)
    decision_det = det(decision_mat)
    if sign(decision_det) != (-1)^(length(S))
        d_nonzero .*= -1
    end

    d[d_nonzero_idxs] = d_nonzero
    return d
end

function good_exit(point, i, conf::Config)
    Vi = get_Vi(conf, point, i)
    xi = point[i]
    zs = abs(Vi) <= conf.ϵ

    return zs || (xi == 0 && Vi < conf.ϵ) || (xi == 1 && Vi > -conf.ϵ), zs
end

function bad_exit(point, i, S, direction)
    tr_coords = vcat(S, [i])
    point_tr = point[tr_coords]
    direction_tr = direction[tr_coords]

    for (j, (pt, dt)) in enumerate(zip(point_tr, direction_tr))
        if (pt == 0 && dt < 0) || (pt == 1 && dt > 0)
            return tr_coords[j]
        end
    end

    return 0
end

function middling_exit(point, i, S, direction, conf::Config)
    new_point = point + conf.γ*direction
    tr_coords = setdiff(1:i-1, S)
    point_tr = point[tr_coords]
    V_tr = get_V_tr(conf, new_point, tr_coords)
    for (j, (pt, Vt)) in enumerate(zip(point_tr, V_tr))
        if (pt == 0 && Vt > 0) || (pt == 1 && Vt < 0)
            return tr_coords[j]
        end
    end

    return 0
end

# sufficient condition for solution to speed up code
# function is_solution(conf::Config, point; α = 0.1)
#     for i in 1:conf.n
#         Vi = get_Vi(conf, point, i)
#         cond1 = Vi == 0
#         cond2 = Vi > 0 && 1 <= α/(conf.n*Vi)+point[i]
#         cond3 = Vi < 0 && 0 <= -α/(conf.n*Vi)-point[i]
#         if !cond1 && !cond2 && !cond3
#             return false
#         end
#     end
#     return true
# end

# check if point satisfies the variational inequality
function is_solution_VI(point, conf::Config; α = 0.1)
    V_x = get_V(conf, point)
    y_worst = ifelse.(V_x .<= 0, 0, 1)
    dot1 = dot(V_x, point)
    dot2 = dot(V_x,y_worst)
    return dot1 - dot2 >= -α
end

"""
run_dynamics(conf::Config)

Executes STON'R dynamics with a given configuration conf.

# Arguments
- `conf::Config` : See `config_FD`, `Config_sym`.

# Output
- `point`: The min-max critical point.
- `pts`: A vector of trajectory points explored during the dynamics.
- `m`: The number of epochs (outer iterations) performed.
- `k`: The number of iterations within each epoch.

# Example

```julia-repl
julia> n = 2
julia> min_coords = [1]
julia> γ = 1e-3
julia> ϵ = 1e-2
julia> f_fd(x) = -2*x[1]*x[2]^2+x[1]^2+x[2]
julia> domain = Hyperrectangle([[-1,1],[-1,1]])
julia> conf = Config_FD(f_fd, n, min_coords, γ, ϵ, domain)
julia> min_max, trajectory, m, k = run_dynamics(conf)
```
"""
function run_dynamics(conf::Config)
    point = fill(0, conf.n)
    i = 1
    S = Vector{Int}()
    pts = Vector{Vector{Float64}}()
    m = 0
    k = 0

    while  i <= conf.n && !is_solution_VI(point, conf)
        k = 0
        println("starting epoch ($i, $S) at $point")

        while true
            direction = compute_direction(point, i, S, conf)
            is_good_exit, is_zs = good_exit(point, i, conf)
            if is_good_exit
                println("good exit at $point")
                if is_zs
                    S = vcat(S,[i])
                end
                i += 1
                break
            end
            
            j = bad_exit(point, i, S, direction)
            if j != 0
                println("bad exit at $point")
                if j == i
                    S = setdiff(S, [i-1])
                    i -= 1
                else
                    S = setdiff(S, [j])
                end
                break
            end

            j = middling_exit(point, i, S, direction, conf)
            if j != 0
                println("middling exit at $point")
                S = vcat(S, [j])
                break
            end
            
            point += conf.γ*direction
            point = P(point)

            k += 1
            push!(pts, point)
        end

        m += 1
    end
    return conf.domain isa Default ? (point, pts, m, k) :
        (H(conf.domain)(point), map(x->H(conf.domain)(x),pts), m, k)
end