# StayOnTheRidge.jl
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://gitlab.fel.cvut.cz/kosohmar/StayOnTheRidge.jl/-/blob/main/LICENSE)

This package implements the STay-ON-the-Ridge algorithm (STON'R), designed to find min-max critical points of possibly nonconvex-nonconcave functions. The concepts of min-max critical points and the STON'R algorithm were proposed in https://proceedings.mlr.press/v195/daskalakis23b.html.

This is a project in progress, developed as part of a Bachelor's thesis at the Faculty of Electrical Engineering, Czech Technical University in Prague.

## Instalation
The package is not registered and this can be installed in the following way

```julia
(@v1.9) pkg> add git@github.com:B0B36JUL-FinalProjects-2023/Projekt_Kosohorska.git
```

## Description
The STON'R algorithm involves the computation of the gradient and hessian of the function. This implementation is able to switch between symbolic computation (using Symbolics.jl) and automatic differentiation (using ForwardDiff.jl).

The algorithm operates on the unit hypercube $\[0,1\]^n$. The domain can be extended to the general convex set $K$ by defining mapping $H$ from $\[0,1\]^n$ to $K$.
The current implementation can operate on $\[0,1\]^n$ or on the general hyperrectangle.

## Example
### Function $2x_1x_2^2 - x_1^2 - x_2$ on the square $\[-1,1\]^2$
To execute dynamics with the `run_dynamics` function, define configuration object first. Two configuration types are `Config_FD` and `Config_sym`. See help for Config_FD and Config_sym for more information.

Execution using ForwardDiff differentiation:
```julia
using StayOnTheRidge

f(x) = 2*x[1]*x[2]^2-x[1]^2-x[2]
n = 2
min_coords = [2]
γ = 1e-3
ϵ = 1e-1
domain = Hyperrectangle([[-1,1],[-1,1]])

config = Config_FD(f, n, min_coords, γ, ϵ, domain)
elapsed = @elapsed min_max, trajectory, m, k = run_dynamics(config)
pretty_print(min_max, elapsed, m, k)
plot_trajectory2D(min_max, trajectory, config.domain)
plot_contour2D(min_max, config.f, config.domain)
plot_surface(min_max, config.f, config.domain)
```

<p align="center">
  <img src="imgs/example5_trajectory.png">
</p>
<p align="center">
  <img src="imgs/example5_contour.png">
</p>

For more examples see examples/examples.jl

