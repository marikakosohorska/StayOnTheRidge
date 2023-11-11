module Stonr

using Symbolics
using LinearAlgebra
using Plots

export execute_stonr, Settings, plot_trajectory2D

struct Settings
    x
    n
    min_coords
    γ
    ϵ
    cube_lb
    cube_rb
    determine_interval
    V
    J
    function Settings(f, x, n, min_coords, γ, ϵ, cube_lb, cube_rb, determine_interval)
        Vs = map(e->prepare_V(e,x,min_coords),f)
        Js = map(e->prepare_J(e,x),Vs)
        return new(x, n, min_coords, γ, ϵ, cube_lb, cube_rb, determine_interval, Vs, Js)
    end
end

function get_V(s::Settings, p)
    return s.V[s.determine_interval(p)]
end

function get_J(s::Settings, p)
    return s.J[s.determine_interval(p)]
end

function prepare_V(f, x, min_coords)
    V = Symbolics.gradient(f, x)
    V[min_coords] .*= -1
    return V
end

function prepare_J(V, x)
    J = Symbolics.jacobian(V, x)
    return J
end

_eval_expr(expr, dict) = Symbolics.value.(Symbolics.substitute(expr, dict))

function compute_direction(point::Vector, i::Integer, S::Vector, s::Settings)
    d = zeros(s.n)

    if isempty(S)
        d[i] = 1
        return d
    end

    d_nonzero_indices = vcat(S, [i])
    #println("d_nonzero_indices $d_nonzero_indices")
    J = get_J(s, point)
    jacobian_Si = J[S, d_nonzero_indices]
    #println("jacobian_Si symbolics: $jacobian_Si")
    jacobian_Si_at_point = _eval_expr(jacobian_Si, Dict(s.x .=> point))
    #println("jacobian_Si_at_point: $jacobian_Si_at_x")
    d_nonzero = nullspace(jacobian_Si_at_point)
    #println("d_nonzero base: $d_nonzero")
    if size(d_nonzero, 2) != 1
        error("Assumption 1 violated at i = $i, S = $S, x = $point, nullspace dimension is $(ndims(d_nonzero)), jacobian_Si_at_point = $jacobian_Si_at_point")
    end

    d_nonzero = vec(d_nonzero)
    decision_matrix = hcat(transpose(jacobian_Si_at_point), d_nonzero)
    #println("decision_matrix $decision_matrix")
    decision_determinant = det(decision_matrix)
    #println("decision_determinant $decision_determinant")
    if sign(decision_determinant) != (-1)^(length(S))
        d_nonzero .*= -1
    end

    d[d_nonzero_indices] = d_nonzero
    #println("Direction is $d")
    return d
end

function good_exit(point::Vector, i::Integer, s::Settings)
    V = get_V(s,point)
    Vi = _eval_expr(V[i],Dict(s.x .=> point))
    xi = point[i]
    zs = abs(Vi) <= s.ϵ

    return zs || (xi == s.cube_lb && Vi < s.ϵ) || (xi == s.cube_rb && Vi > -s.ϵ), zs
end

function bad_exit(point::Vector, i::Integer, S::Vector, direction::Vector, s::Settings)
    triggering_coords = vcat(S, [i])
    point_triggering = point[triggering_coords]
    direction_triggering = direction[triggering_coords]

    for (j, (pt, dt)) in enumerate(zip(point_triggering, direction_triggering))
        if (pt == s.cube_lb && dt < 0) || (pt == s.cube_rb && dt > 0)
            return triggering_coords[j]
        end
    end

    return 0
end

function middling_exit(point::Vector, i::Integer, S::Vector, direction::Vector, s::Settings)
    new_point = point + s.γ*direction
    triggering_coords = setdiff(1:i-1, S)
    point_triggering = point[triggering_coords]
    V = get_V(s,point)
    V_triggering = _eval_expr(V[triggering_coords], Dict(s.x .=> new_point))
    for (j, (pt, Vt)) in enumerate(zip(point_triggering, V_triggering))
        if (pt == s.cube_lb && Vt > 0) || (pt == s.cube_rb && Vt < 0)
            return triggering_coords[j]
        end
    end

    return 0
end

P(x; x_min = 0, x_max = 1) = min.(max.(x, x_min), x_max)

function execute_stonr(s::Settings)
    # init values
    point = fill(s.cube_lb,s.n)
    i = 1
    S = Vector{Int}()
    pts = []

    while i <= s.n
        println("starting epoch ($i, $S) at $point")
        while true
            is_good_exit, is_zs = good_exit(P(point; x_min = s.cube_lb, x_max = s.cube_rb), i, s)
                if is_good_exit
                    println("good exit at $point")
                    if is_zs
                        S = vcat(S,[i])
                    end
                    i += 1
                    break
                end
            direction = compute_direction(point, i, S, s)
            j = bad_exit(P(point; x_min = s.cube_lb, x_max = s.cube_rb), i, S, direction,s)
            if j != 0
                println("bad exit at $point")
                if j == i
                    i -= 1
                    S = setdiff(S, [i-1])
                else
                    S = setdiff(S, [j])
                end
                break
            end
            j = middling_exit(P(point; x_min = s.cube_lb, x_max = s.cube_rb), i, S, direction,s)
            if j != 0
                println("middling exit at $point")
                S = vcat(S,[j])
                break
            end
            point = point + s.γ*direction
            direction = compute_direction(point, i, S,s)
            push!(pts, point)
        end
        point = P(point; x_min = s.cube_lb, x_max = s.cube_rb)
    end

    return point, pts
end

function plot_trajectory2D(min_max::Vector, trajectory, s::Settings)
    x1_coords = [pt[1] for pt in trajectory]
    x2_coords = [pt[2] for pt in trajectory]
    scatter(x1_coords, x2_coords; color=:blue, legend = false, xlims = (s.cube_lb,s.cube_rb), ylims = (s.cube_lb,s.cube_rb), markersize = 1, markershape=:auto)
    scatter!([min_max[1]], [min_max[2]], markershape=:star, markersize=10)
end

end
