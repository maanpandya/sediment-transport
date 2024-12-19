using HomotopyContinuation
using LinearAlgebra
using DynamicPolynomials

function solve_polynomial_system(input_alpha, input_beta, input_gamma, input_delta, input_omega, input_funcs, returnnonsingular=false)
    n = length(input_alpha) * 2  # number of unknown variables

    # Declare variables
    @polyvar u[1:n÷2] v[1:n÷2]

    # Assign parameters
    α = input_alpha
    β = input_beta
    γ = input_gamma
    δ = input_delta
    ω = input_omega

    # Define the system of equations
    f = input_funcs

    # Define the system as a polynomial system
    system = HomotopyContinuation.System(f)

    # Solve the system
    result = HomotopyContinuation.solve(system)

    # Collect real solutions
    real_solutions = []
    for sol in result
        if HomotopyContinuation.is_real(sol)
            push!(real_solutions, sol)
        end
    end

    # Return nonsingular solutions if requested
    if returnnonsingular
        nonsingular_solutions = []
        for sol in real_solutions
            if HomotopyContinuation.is_nonsingular(sol)
                push!(nonsingular_solutions, sol)
            end
        end
        return real_solutions, nonsingular_solutions
    end

    return real_solutions
end

n = 4 # number of unknown variables

# declare the variables 
@polyvar u[1:n÷2] v[1:n÷2]
@polyvar f[1:n]

# declare the parameters
@polyvar α[1:n÷2]
@polyvar β[1:n÷2]
@polyvar γ[1:n÷2]
@polyvar δ[1:n÷2]

# this input funcs will take Hew's function as an input containing all the polynomial equations
input_alpha = [1.0, 1.1]
input_beta = [0.04, 0.03]
input_gamma = [1.0, 1.1]
input_delta = [0.1, 0.2]
input_omega = [1.5]

α = input_alpha
β = input_beta
γ = input_gamma
δ = input_delta
ω = input_omega

input_function = [
    (ω[1]^2 - α[1])*u[1] - δ[1]*ω[1]*v[1] + γ[1] - (3/4)*β[1]*(u[1]^3 + u[1]*v[1]^2),
    (ω[1]^2 - α[2])*u[2] - δ[2]*ω[1]*v[2] + γ[2] - (3/4)*β[2]*(u[2]^3 + u[2]*v[2]^2),
    ω[1]*δ[1]*u[1] + (ω[1]^2 - α[1])*v[1] - (3/4)*β[1]*(v[1]^3 + u[1]^2*v[1]),
    ω[1]*δ[2]*u[2] + (ω[1]^2 - α[2])*v[2] - (3/4)*β[2]*(v[2]^3 + u[2]^2*v[2])
]

solve_polynomial_system(input_alpha, input_beta, input_gamma, input_delta, input_omega, input_function)