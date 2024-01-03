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
    return Config_sym(f, x, n, min_coords, γ, ϵ, [[0,1],[0,1]])
end

function conf2()
    n = 4
    min_coords = [1,2]
    γ = 1e-3
    ϵ = 1e-2
    x = [Symbolics.variable(:x,i) for i in 1:n]
    f = (x[1]+x[2]+x[3]+x[4]+1)^2-4*(x[1]*x[2]+x[2]*x[3]+x[3]*x[4]+x[1])
    return Config_sym(f, x, n, min_coords, γ, ϵ, [[0,1],[0,1],[0,1],[0,1]])
end

function conf3()
    n = 6
    min_coords = [1,2,3]
    γ = 1e-4
    ϵ = 1e-1
    f(x) = sum(x[i]+x[3+i] for i in 1:3) - prod((x[i] - x[3+i]) for i in 1:3)
    return Config_FD(f, n, min_coords, γ, ϵ, [fill([-1,1], 6)...])
end

@testset "StayOnTheRidge.jl" begin
   c1 = conf1()
   c2 = conf2()
   c3 = conf3()
   @testset "Direction" begin
        @test StayOnTheRidge.compute_direction([0,0], 1, [], c1) == [1,0]
        @test StayOnTheRidge.compute_direction([0,0.5], 1, [], c1) == [1,0]
        @test StayOnTheRidge.compute_direction([1,0], 2, [], c1) == [0,1]
        @test StayOnTheRidge.compute_direction([1,0.5], 2, [], c1) == [0,1]
        @test StayOnTheRidge.compute_direction([1,0.5], 2, [1], c1) == [-1,0]
        @test_throws ErrorException StayOnTheRidge.compute_direction([1,0,0,0], 3, [1,2], c2)
   end
   @testset "Exit points" begin
    is_good_exit,zs = StayOnTheRidge.good_exit([1,0],1,c1)
    @test is_good_exit && !zs
    j = StayOnTheRidge.middling_exit([1,0.5],2,[],[0,1],c1)
    @test j == 1
    is_good_exit,zs = StayOnTheRidge.good_exit([0.5,0.5],2,c1)
    @test is_good_exit && zs
   end
   @testset "VI solution" begin
    @test StayOnTheRidge.is_solution_VI([-1,-1,1,1,1,1],c3)
    @test StayOnTheRidge.is_solution_VI([-1,1,-1,1,1,1],c3)
    @test StayOnTheRidge.is_solution_VI([1,-1,-1,1,1,1],c3)
    @test !StayOnTheRidge.is_solution_VI([0,0,0,0,0,0],c3)
   end
end
