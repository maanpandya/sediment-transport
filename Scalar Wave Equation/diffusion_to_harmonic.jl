using Symbolics, HarmonicBalance, OrderedCollections

# Define variables
@variables t
@variables (u(t))[1:3]  # Variables for interior nodes u1, u2, u3

# Spatial step size
Δx = 0.25

# Equations of motion
eqs = [
    d(u[1], t) ~ (0 - 2u[1] + u[2]) / Δx^2,
    d(u[2], t) ~ (u[1] - 2u[2] + u[3]) / Δx^2,
    d(u[3], t) ~ (u[2] - 2u[3] + 0) / Δx^2
]

# Flatten u to a vector of variables
vars = collect(u)  # Convert Symbolics.Arr to Vector{Num}

# Create the DifferentialEquation object
diff_eq = DifferentialEquation(eqs, vars)

# Assign harmonics to variables
for i in vars
    add_harmonic!(diff_eq, i, [1])  # Add the fundamental harmonic
end

# Display the differential equation object
println(diff_eq)
harmonic_eq = get_harmonic_equations(diff_eq)
println(harmonic_eq)
