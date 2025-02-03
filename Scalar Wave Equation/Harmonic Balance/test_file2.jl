test_file2.jl
using TestSetExtensions

@testset "Test Set Extensions for test_file2" begin
    @testset "Basic Tests" begin
        @test 1 + 1 == 2
        @test 2 * 2 == 4
    end

    @testset "Advanced Tests" begin
        @test 3 - 1 == 2
        @test 4 / 2 == 2
    end
end