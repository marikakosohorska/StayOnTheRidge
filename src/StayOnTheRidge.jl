module StayOnTheRidge

using LinearAlgebra
using Plots
using ForwardDiff
using Symbolics

export run_dynamics, Config_sym, Config_FD, plot_trajectory2D, pretty_print, format_elapsed

abstract type Config end

struct Config_sym <: Config
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

    function Config_sym(f_expr, x, n, min_coords, γ, ϵ; H = nothing)
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

struct Config_FD <: Config
    f
    n
    min_coords
    γ
    ϵ
    H
    f_obj

    function Config_FD(f, n, min_coords, γ, ϵ; H = nothing)
        f_obj = isnothing(H) ? f : x -> f(H(x))
        return new(f, n, min_coords,  γ, ϵ, H, f_obj)
    end
end

eval_expr(expr, dict) = Symbolics.value.(Symbolics.substitute(expr, dict))

function prepare_V(f, x, min_coords)
    V = Symbolics.gradient(f, x)
    V[min_coords] .*= -1
    return V
end

function prepare_V(f, x, min_coords, H)
    f_at_H = eval_expr(f, Dict(x .=> H(x)))
    V = Symbolics.gradient(f_at_H, x)
    V[min_coords] .*= -1
    return V
end

P(x; x_min = 0, x_max = 1) = min.(max.(x, x_min), x_max)

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

function get_V_tr(conf::Config_sym, new_point, tr_coords)
    V_tr = [conf.V[tr](new_point) for tr in tr_coords]
    return V_tr
end

function get_V_tr(conf::Config_FD, new_point, tr_coords)
    V_tr = get_V(conf, new_point)[tr_coords]
    return V_tr
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

function is_solution(conf::Config, point) # sufficient condition for solution to speed up code
    α = 0.1
    for i in 1:conf.n
        Vi = get_Vi(conf, point, i)
        cond1 = Vi == 0
        cond2 = Vi > 0 && 1 <= α/(conf.n*Vi)+point[i]
        cond3 = Vi < 0 && 0 <= -α/(conf.n*Vi)-point[i]
        if !cond1 && !cond2 && !cond3
            return false
        end
    end
    return true
end

function run_dynamics(conf::Config)
    point = fill(0, conf.n)
    i = 1
    S = Vector{Int}()
    pts = Vector{Vector{Float64}}()
    m = 0
    k = 0

    while i <= conf.n && !is_solution(conf, point)
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

    return point, pts, m, k
end

function plot_trajectory2D(min_max, trajectory, min_bound, max_bound)
    x1_coords = [pt[1] for pt in trajectory]
    x2_coords = [pt[2] for pt in trajectory]
    scatter(x1_coords, x2_coords; color=:blue, legend = false, markersize = 1,
        markershape=:auto, xlims = [min_bound-0.1, max_bound+0.1], ylims = [min_bound-0.1, max_bound+0.1])
    scatter!([min_max[1]], [min_max[2]], markershape=:star, markersize=10)
end

function pretty_print(point, time, m, k)
    str_point = "($(join(round.(point, digits=4), ", ")))"
    str_time = format_elapsed(time)
    display("Min max critical point: $str_point, elapsed time: $str_time, outer cycles: $m, inner cycles: $k")
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

end


