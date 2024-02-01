using Plots

function get_bounds(domain::Hyperrectangle)
    p11 = domain.sides[1][1]
    p12 = domain.sides[1][2]
    p21 = domain.sides[2][1]
    p22 = domain.sides[2][2]
    x1_min_bound = p11 - abs(p11-p12)/20
    x1_max_bound = p12 + abs(p11-p12)/20
    x2_min_bound = p21 - abs(p21-p22)/20
    x2_max_bound = p22 + abs(p21-p22)/20
    return [x1_min_bound, x1_max_bound, x2_min_bound, x2_max_bound]
end

function get_bounds(domain::Default)
    return [-1/20,21/20,-1/20,21/20]
end

function get_ranges(domain::Hyperrectangle)
    return [range(domain.sides[1][1], domain.sides[1][2], length=300),
        range(domain.sides[2][1], domain.sides[2][2], length=300)]
end

function get_ranges(domain::Default)
    return [range(0, 1, length=300),
    range(0, 1, length=300)]
end

function plot_trajectory2D(min_max, trajectory, domain::Domain)
    x1_coords = [pt[1] for pt in trajectory]
    x2_coords = [pt[2] for pt in trajectory]
    bounds = get_bounds(domain)
    scatter(x1_coords, x2_coords; color=:blue, legend = false, markersize = 1,
        markershape=:auto, xlims = [bounds[1], bounds[2]], ylims = [bounds[3], bounds[4]])
    scatter!([min_max[1]], [min_max[2]], markershape=:star, markersize=10, color=:blue)
    title!("Trajectory")
end

function plot_contour2D(min_max, f, domain::Domain)
    f_plot(x,y) = f([x,y])
    rs = get_ranges(domain)
    bounds = get_bounds(domain)
    contourf(rs[1], rs[2], f_plot, levels=50, xlims = [bounds[1], bounds[2]], ylims = [bounds[3], bounds[4]])
    scatter!([min_max[1]], [min_max[2]], markershape=:star, markersize=10, legend=false, color=:blue)
    title!("Contour plot of the function")
end

function plot_surface(min_max, f, domain::Domain)
    f_plot(x,y) = f([x,y])
    rs = get_ranges(domain)
    zs = [f_plot(x,y) for y in rs[2], x in rs[1]]
    surface(rs[1], rs[2], zs, legend=false)
    scatter!([min_max[1]], [min_max[2]], [f(min_max)], markershape=:star, markersize=10, legend=false, color=:blue)
    title!("3D plot of the function")
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