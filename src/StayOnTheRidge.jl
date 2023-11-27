module StayOnTheRidge

using Symbolics
using LinearAlgebra
using Plots

export run_dynamics, Settings, plot_trajectory2D

struct Settings
    f
    x
    n
    min_coords
    γ
    ϵ
    V
    J
    grad_f
    hessian_f
    H

    function Settings(f_expr, x, n, min_coords, γ, ϵ; H = nothing)
        f = eval(build_function(f_expr, x))
        V_expr = isnothing(H) ? prepare_V(f_expr, x, min_coords) : prepare_V(f_expr, x, min_coords, H)
        V = [eval(build_function(V_expr[i], x)) for i in 1:n]
        J_expr = Symbolics.jacobian(V_expr, x)
        J = [eval(build_function(J_expr[i,j], x)) for i in 1:n, j in 1:n]
        grad_f_expr = Symbolics.gradient(f_expr, x)
        grad_f = [eval(build_function(grad_f_expr[i], x)) for i in 1:n]
        hessian_f_expr = Symbolics.jacobian(grad_f_expr, x)
        hessian_f = [eval(build_function(hessian_f_expr[i,j], x)) for i in 1:n, j in 1:n]
        return new(f, x, n, min_coords, γ, ϵ, V, J, grad_f, hessian_f, H)
    end
end

_eval_expr(expr, dict) = Symbolics.value.(Symbolics.substitute(expr, dict))

function prepare_V(f, x, min_coords)
    V = Symbolics.gradient(f, x)
    V[min_coords] .*= -1
    return V
end

function prepare_V(f, x, min_coords, H)
    f_at_H = _eval_expr(f, Dict(x .=> H(x)))
    V = Symbolics.gradient(f_at_H, x)
    V[min_coords] .*= -1
    return V
end

P(x; x_min = 0, x_max = 1) = min.(max.(x, x_min), x_max)

function compute_direction(point::Vector, i::Integer, S::Vector, s::Settings)
    d = zeros(s.n)

    if isempty(S)
        d[i] = 1
        return d
    end

    d_nonzero_idxs = vcat(S, [i])
    jacobian_Si = [s.J[i,j](point) for i in S, j in d_nonzero_idxs]
    d_nonzero = nullspace(jacobian_Si)
    if size(d_nonzero, 2) != 1
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
 
function good_exit(point::Vector, i::Integer, s::Settings)
    Vi = s.V[i](point)
    xi = point[i]
    zs = abs(Vi) <= s.ϵ

    return zs || (xi == 0 && Vi < s.ϵ) || (xi == 1 && Vi > -s.ϵ), zs
end

function bad_exit(point::Vector, i::Integer, S::Vector, direction::Vector)
    triggering_coords = vcat(S, [i])
    point_triggering = point[triggering_coords]
    direction_triggering = direction[triggering_coords]

    for (j, (pt, dt)) in enumerate(zip(point_triggering, direction_triggering))
        if (pt == 0 && dt < 0) || (pt == 1 && dt > 0)
            return triggering_coords[j]
        end
    end

    return 0
end

function middling_exit(point::Vector, i::Integer, S::Vector, direction::Vector, s::Settings)
    new_point = point + s.γ*direction
    triggering_coords = setdiff(1:i-1, S)
    point_triggering = point[triggering_coords]
    V_triggering = [s.V[tr](new_point) for tr in triggering_coords]
    for (j, (pt, Vt)) in enumerate(zip(point_triggering, V_triggering))
        if (pt == 0 && Vt > 0) || (pt == 1 && Vt < 0)
            return triggering_coords[j]
        end
    end

    return 0
end

function run_dynamics(s::Settings)
    point = fill(0, s.n)
    i = 1
    S = Vector{Int}()
    pts = Vector{Vector{Float64}}()

    m = 0
    k = 0
    while i <= s.n
        m += 1
        k = 0
        println("starting epoch ($i, $S) at $point")

        while true
            k += 1
            proj = P(point)
            direction_proj = compute_direction(proj, i, S, s)

            is_good_exit, is_zs = good_exit(proj, i, s)
            if is_good_exit
                println("good exit at $point")
                if is_zs
                    S = vcat(S,[i])
                end
                i += 1
                break
            end

            j = bad_exit(proj, i, S, direction_proj)
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

            j = middling_exit(proj, i, S, direction_proj, s)
            if j != 0
                println("middling exit at $point")
                S = vcat(S, [j])
                break
            end
            
            direction = compute_direction(point, i, S, s)
            point += s.γ*direction
            push!(pts, point)
        end

        point = P(point)
    end

    return point, pts, m, k
end

function plot_trajectory2D(min_max::Vector, trajectory, min_bound, max_bound)
    x1_coords = [pt[1] for pt in trajectory]
    x2_coords = [pt[2] for pt in trajectory]
    scatter(x1_coords, x2_coords; color=:blue, legend = false, markersize = 1,
        markershape=:auto, xlims = [min_bound-0.1, max_bound+0.1], ylims = [min_bound-0.1, max_bound+0.1])
    scatter!([min_max[1]], [min_max[2]], markershape=:star, markersize=10)
end

end
