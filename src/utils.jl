using Plots

function plot_trajectory2D(min_max, trajectory, min_bound, max_bound)
    x1_coords = [pt[1] for pt in trajectory]
    x2_coords = [pt[2] for pt in trajectory]
    scatter(x1_coords, x2_coords; color=:blue, legend = false, markersize = 1,
        markershape=:auto, xlims = [min_bound-0.1, max_bound+0.1], ylims = [min_bound-0.1, max_bound+0.1])
    scatter!([min_max[1]], [min_max[2]], markershape=:star, markersize=10, color=:blue)
    title!("Trajectory")
end

function plot_contour2D(min_max, f, min_bound, max_bound)
    f_plot(x,y) = f([x,y])
    x1 = range(min_bound, max_bound, length=300)
    x2 = range(min_bound, max_bound, length=300)
    contourf(x1, x2, f_plot, levels=50)
    scatter!([min_max[1]], [min_max[2]], markershape=:star, markersize=10, legend=false, color=:blue)
    title!("Contour plot of the function")
end

function plot_surface(min_max, f, min_bound, max_bound)
    f_plot(x,y) = f([x,y])
    xs = ys = range(min_bound, max_bound, length=100)
    zs = [f_plot(x,y) for y in ys, x in xs]
    surface(xs, ys, zs, legend=false)
    scatter!([min_max[1]], [min_max[2]], [f(min_max)], markershape=:star, markersize=10, legend=false, color=:blue)
    title!("3D plot of the function")
end

# function plot_function(min_max, min_bound, max_bound, f)
#     p1 = plot_contour2D(min_max, min_bound, max_bound, f)
#     p2 = plot_surface(min_max, min_bound, max_bound, f)
#     l = @layout [a{0.5w} ; b]
#     plot(p1, p2, layout=l)
# end

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