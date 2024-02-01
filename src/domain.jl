abstract type Domain end

"""
Default domain [0,1]ⁿ.
"""
struct Default <: Domain end # unit hypercube [0,1]ⁿ

"""
Hyperrectangle

# Arguments
- `sides::Vector{Vector{Real}}` : Define the hyperrectangle, e.g. sides = [[-1,1],[0,1]] define the rectangle [-1,1]×[0,1].
"""
struct Hyperrectangle <: Domain
    sides::Vector{Vector{T}} where T <: Real
    function Hyperrectangle(sides::Vector{Vector{T}} where T <: Real)
        @assert all(subvec -> length(subvec) >= 2 && subvec[1] < subvec[2], sides)
        return new(sides)
    end
end

# general hyperrectangle mapping
function H(domain::Hyperrectangle)
    sides = domain.sides
    return x -> getindex.(sides,1) .+ (getindex.(sides,2) .- getindex.(sides,1)) .* x
end