using Test
using TestSetExtensions

@testset "Harmonic Balance Tests" begin
    @testset "Test Case 1" begin
        @test 1 + 1 == 2
    end

    @testset "Test Case 2" begin
        @test 2 * 2 == 4
    end
end

@testset "Additional Tests" begin
    @testset "Test Case 3" begin
        @test 3 - 1 == 2
    end

    @testset "Test Case 4" begin
        @test 5 / 5 == 1
    end
end