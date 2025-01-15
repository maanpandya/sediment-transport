using HomotopyContinuation
using LinearAlgebra
using DynamicPolynomials
using Symbolics
using SymbolicUtils

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
@variables u[1:n÷2] v[1:n÷2]
@variables f[1:n]

# declare the parameters
@variables α[1:n÷2]
@variables β[1:n÷2]
@variables γ[1:n÷2]
@variables δ[1:n÷2]

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


println(typeof(input_function))
new_eq = input_function[1]

println(new_eq)
println(typeof(new_eq))

@polyvar u_poly[1:1] v_poly[1:1]

# Parameter values
ω_val = [2.0]
α_val = [1.0]
δ_val = [0.1]
γ_val = [1.0]
β_val = [0.04]

# Create polynomial equation by substituting numeric values and polynomial variables
poly_eq = Symbolics.substitute(new_eq, 
    Dict(
        # ω[1] => ω_val[1],
        # α[1] => α_val[1],
        # δ[1] => δ_val[1],
        # γ[1] => γ_val[1],
        # β[1] => β_val[1],
        u[1] => u_poly[1],
        v[1] => v_poly[1]
    )
)

println(typeof(poly_eq))

# result = solve(system)
# println(result)

# solve_polynomial_system(input_alpha, input_beta, input_gamma, input_delta, input_omega, input_function)