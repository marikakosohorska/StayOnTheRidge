abstract type Config end

"""
Config_FD

Configuration structure to store data for STON'R execution.
Using ForwardDiff.jl for differentiation.

# Arguments
- `f::Function` : Function on domain [0,1]ⁿ.
- `n::Int` : Number of variables.
- `min_coords::Vector{Int}` : Vector of minimizing coordinates.
- `γ::Real` : Step size.
- `ϵ::Real` : Precision parameter of min-max critical point.
- `domain::Domain` : Specifies the domain of the function (see Default or Hyperrectangle).
# Example
```julia-repl
julia> f_fd(x) = -2*x[1]*x[2]^2+x[1]^2+x[2]
julia> n = 2
julia> min_coords = [1]
julia> γ = 1e-2
julia> ϵ = 1e-2
julia> domain = Default()
julia> conf = Config_FD(f_fd, n, min_coords, γ, ϵ, domain);
```
"""
struct Config_FD <: Config
    f::Function
    f_obj::Function
    n::Int
    min_coords::Vector{Int}
    γ::Real
    ϵ::Real
    domain::Domain

    function Config_FD(f::Function, n::Int, min_coords::Vector{Int},
        γ::Real, ϵ::Real, domain::Domain)
        @assert n > 0
        @assert all(x -> x > 0, min_coords)
        @assert γ > 0
        @assert ϵ > 0
        @assert length(min_coords) <= n
        if domain isa Default # default [0,1]ⁿ
            f_obj = f
        else
            inner = H(domain) # closure to gain performance
            f_obj = x -> f(inner(x))
        end
        return new(f, f_obj, n, min_coords, γ, ϵ, domain)
    end
end

"""
Config_sym

Configuration structure to store data for STON'R execution.
Using Symbolics.jl for differentiation.

# Arguments
- `f_expr::Num` : Symbolic function on domain [0,1]ⁿ.
- `x::Num` : Vector of symbolic variables.
- `n::Int` : Number of variables.
- `min_coords::Vector{Int}` : Vector of minimizing coordinates.
- `γ::Real` : Step size.
- `ϵ::Real` : Precision parameter of min-max critical point.
- `domain::Domain` : Specifies the domain of the function (see Default or Hyperrectangle).
# Example
```julia-repl
julia> n = 2
julia> x = [Symbolics.variable(:x,i) for i in 1:n]
julia> f = (x[1]-1/2)*(x[2]-1/2)
julia> min_coords = [1]
julia> γ = 1e-2
julia> ϵ = 1e-2
julia> domain = Default()
julia> conf = Config_sym(f, x, n, min_coords, γ, ϵ, domain);
```
"""
struct Config_sym <: Config
    f::Function
    V::Vector{Function} # modified gradient of f
    J::Matrix{Function} # modified Jacobian matrix of V
    x::Vector{Num}
    n::Int
    min_coords::Vector{Int}
    γ::Real
    ϵ::Real
    domain::Domain

    # for testing
    grad_f::Vector{Function}
    hess_mat_f::Matrix{Function}

    function Config_sym(f_expr::Num, x::Vector{Num}, n::Int, min_coords::Vector{Int},
        γ::Real, ϵ::Real, domain::Domain)
        @assert n > 0
        @assert all(x -> x > 0, min_coords)
        @assert γ > 0
        @assert ϵ > 0
        @assert length(x) == n
        @assert length(min_coords) <= n
        f = eval(build_function(f_expr, x))
        V_expr = prepare_V(f_expr, x, min_coords, domain)
        V = [eval(build_function(V_expr[i], x)) for i in 1:n]
        J_expr = Symbolics.jacobian(V_expr, x)
        J = [eval(build_function(J_expr[i,j], x)) for i in 1:n, j in 1:n]
        grad_f_expr = Symbolics.gradient(f_expr, x)
        grad_f = [eval(build_function(grad_f_expr[i], x)) for i in 1:n]
        hess_mat_f_expr = Symbolics.jacobian(grad_f_expr, x)
        hess_mat_f = [eval(build_function(hess_mat_f_expr[i,j], x)) for i in 1:n, j in 1:n]
        return new(f, V, J, x, n, min_coords, γ, ϵ, domain, grad_f, hess_mat_f)
    end
end

eval_expr(expr, dict) = Symbolics.value.(Symbolics.substitute(expr, dict))

function prepare_V(f, x, min_coords, domain)
    if domain isa Default
        f_obj = f
    else
        f_obj = eval_expr(f, Dict(x .=> H(domain)(x)))
    end
    V = Symbolics.gradient(f_obj, x)
    V[min_coords] .*= -1
    return V
end