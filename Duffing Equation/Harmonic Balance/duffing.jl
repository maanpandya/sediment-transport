# Import necessary libraries
using HarmonicBalance
using Plots

# Define the Duffing equation coefficients
α = 1.0       # Linear stiffness
β = 5.0       # Nonlinear stiffness
δ = 0.02      # Damping coefficient
γ = 8.0       # Driving force amplitude
ω = 1.0       # Driving force frequency

# Function to define the Duffing equation
function duffing_eqn(t, x)
    dxdt = zeros(2)
    dxdt[1] = x[2]
    dxdt[2] = -δ * x[2] - α * x[1] - β * x[1]^3 + γ * cos(ω * t)
    return dxdt
end

# Number of harmonics
n_harmonics = 5  # You can change this value for more harmonics

# Set the initial conditions for displacement and velocity
x_init = 1.0  # Initial displacement
v_init = 0.0  # Initial velocity

# Function to convert initial conditions into Fourier coefficients for Harmonic Balance
function set_initial_conditions(x_init, v_init, n_harmonics)
    x0 = zeros(2 * n_harmonics)  # Initial guess vector
    # Set initial displacement (for the first harmonic, corresponding to the DC component)
    x0[1] = x_init
    # Set initial velocity (for the first harmonic derivative)
    x0[2] = v_init
    return x0
end

# Create the initial guess based on specified initial conditions
x0 = set_initial_conditions(x_init, v_init, n_harmonics)

# Set up the harmonic balance solver
hb_solver = HarmonicBalance(n_harmonics)

# Use the Harmonic Balance method to solve the Duffing equation
result = solve(hb_solver, duffing_eqn, x0)

# Extract the time and solution from the result
t_vals = LinRange(0, 2 * π / ω, 1000)  # Time values for plotting
x_vals = real(time_series(result, t_vals))

# Plot the results
plot(t_vals, x_vals[1, :], xlabel = "Time", ylabel = "Displacement", label = "Duffing Oscillator")
title!("Duffing Equation Solution with Custom Initial Conditions")