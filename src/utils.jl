using Plots

function plot_trajectory2D(min_max, trajectory, domain::Hyperrectangle)
    x1_coords = [pt[1] for pt in trajectory]
    x2_coords = [pt[2] for pt in trajectory]
    p11 = domain.sides[1][1]
    p12 = domain.sides[1][2]
    p21 = domain.sides[2][1]
    p22 = domain.sides[2][2]
    x1_min_bound = p11 - abs(p11-p12)/20
    x1_max_bound = p12 + abs(p11-p12)/20
    x2_min_bound = p21 - abs(p21-p22)/20
    x2_max_bound = p22 + abs(p21-p22)/20
    scatter(x1_coords, x2_coords; color=:blue, aspect_ratio=:equal, legend = false, markersize = 1,
        markershape=:auto, xlims = [x1_min_bound, x1_max_bound], ylims = [x2_min_bound, x2_max_bound])
    scatter!([min_max[1]], [min_max[2]], markershape=:star, markersize=10, color=:blue)
    title!("Trajectory")
end

function plot_contour2D(min_max, f, domain::Hyperrectangle)
    f_plot(x,y) = f([x,y])
    p11 = domain.sides[1][1]
    p12 = domain.sides[1][2]
    p21 = domain.sides[2][1]
    p22 = domain.sides[2][2]
    x1 = range(p11, p12, length=300)
    x2 = range(p21, p22, length=300)
    x1_min_bound = p11 - abs(p11-p12)/20
    x1_max_bound = p12 + abs(p11-p12)/20
    x2_min_bound = p21 - abs(p21-p22)/20
    x2_max_bound = p22 + abs(p21-p22)/20
    contourf(x1, x2, f_plot, levels=50, aspect_ratio=:equal, xlims = [x1_min_bound, x1_max_bound], ylims = [x2_min_bound, x2_max_bound])
    scatter!([min_max[1]], [min_max[2]], markershape=:star, markersize=10, legend=false, color=:blue)
    title!("Contour plot of the function")
end

function plot_surface(min_max, f, domain::Hyperrectangle)
    f_plot(x,y) = f([x,y])
    p11 = domain.sides[1][1]
    p12 = domain.sides[1][2]
    p21 = domain.sides[2][1]
    p22 = domain.sides[2][2]
    xs = range(p11, p12, length=100)
    ys = range(p21, p22, length=100)
    zs = [f_plot(x,y) for y in ys, x in xs]
    surface(xs, ys, zs, legend=false)
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