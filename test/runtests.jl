using StayOnTheRidge
using Test
using Symbolics

function conf1()
    n = 2
    min_coords = [1]
    γ = 1e-2
    ϵ = 1e-2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    f = (x[1]-1/2)*(x[2]-1/2)
    Config_sym(f, x, n, min_coords, γ, ϵ)
end

function conf2()
    n = 4
    min_coords = [1,2]
    γ = 1e-3
    ϵ = 1e-2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    f = (x[1]+x[2]+x[3]+x[4]+1)^2-4*(x[1]*x[2]+x[2]*x[3]+x[3]*x[4]+x[1])
    Config_sym(f, x, n, min_coords, γ, ϵ)
end

@testset "StayOnTheRidge.jl" begin
   c1 = conf1()
   c2 = conf2()
   @testset "Direction" begin
        @test StayOnTheRidge.compute_direction([0,0], 1, [], c1) == [1,0]
        @test StayOnTheRidge.compute_direction([0,0.5], 1, [], c1) == [1,0]
        @test StayOnTheRidge.compute_direction([1,0], 2, [], c1) == [0,1]
        @test StayOnTheRidge.compute_direction([1,0.5], 2, [], c1) == [0,1]
        @test StayOnTheRidge.compute_direction([1,0.5], 2, [1], c1) == [-1,0]
        @test_throws ErrorException StayOnTheRidge.compute_direction([1,0,0,0], 3, [1,2], c2)
   end
end
