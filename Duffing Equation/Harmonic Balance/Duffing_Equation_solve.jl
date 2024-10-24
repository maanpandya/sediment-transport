using LinearAlgebra
using DifferentialEquations
using FFTW
using NLsolve
using BenchmarkTools
using Plots 

# Set initial conditions for position (u0) and velocity (v0)
initial_position = 3.0  # Initial position
initial_velocity = 0.0  # Initial velocity

# Combine into a vector for initial conditions
u0 = [initial_position, initial_velocity]

# Set end time and time interval 
Tend = 20.0  # Total simulation time
tspan = (0.0, Tend)           

# Constants
m = 1.0      # Mass
gamma = 0.1  # Damping coefficient
k1 = 1.0     # Linear stiffness
k3 = 1.0     # Nonlinear stiffness
F0 = 1.0     # Forcing amplitude
omega = 1.0  # Driving frequency

# Define the right-hand side (RHS) of the Duffing equation
function rhs!(ddu, u, p, t)
    F_drive = F0 * sin(omega * t)
    ddu .= zeros(2)  # Ensure ddu is a 2-element vector
    ddu[1] = u[2]  # First derivative (velocity)
    ddu[2] = (-k1*u[1] - gamma*u[2] - k3*u[1]^3 + F_drive) / m  # Second derivative (acceleration)
end

# Define the second-order ODE problem
prob = SecondOrderODEProblem(rhs!, u0, tspan)

# Solve the ODE problem using an adaptive solver
sol = solve(prob, AutoVern7(Rodas5()), reltol=1e-8, abstol=1e-8)

# Set sampling points for the FFT
Nsamples = 500
Tstart = 1.0
dt = (Tend - Tstart) / Nsamples
tsampled = Vector(Tstart:dt:Tend)

# Sample the velocity and position data from the solution
vsampled = [sol(tk, idxs=1) for tk in tsampled]
usampled = [sol(tk, idxs=2) for tk in tsampled]

# Perform FFT on the sampled position data
uf = fft(usampled)

# Set the frequency axis for plotting
fmax = 1 / (2.0 * dt)
fstep = 2 * fmax / Nsamples
fvec = Vector(0:fstep:fmax)

# Plot computed velocity and position over time
p1 = plot(sol, idxs=1, line=:blue, title="Velocity starts at v0")
p2 = plot(sol, idxs=2, line=:red, title="Position starts at u0")
display(plot(p1, p2, layout=(1, 2)))  # Explicitly display the time-domain plots

# Plot the absolute value and phase of FFT samples
p3 = bar(fvec, 2.0 / Nsamples * log.(abs.(uf[1:length(fvec)])), xlabel="Frequency", ylabel="Magnitude (log)", label="Norm")
p4 = plot(fvec, angle.(uf[1:length(fvec)]), xlabel="Frequency", ylabel="Phase", label="Phase")
display(plot(p3, p4, layout=(2, 1)))  # Explicitly display the FFT plots

#----------------------------------------------------------------------------------------------------------------------------------#

# Define the system of equations for harmonic balance
function balance_eqs!(F, x)
    A, B = x[1], x[2]
    F[1] = -m*A*omega^2 + gamma*B*omega + k1*A + (3/4)*k3*A^3 + (3/2)*k3*A*B^2
    F[2] = -m*B*omega^2 - gamma*A*omega + k1*B + (3/4)*k3*B^3 + (3/2)*k3*A^2*B - F0
end

# Initial guess for A and B
x0 = [1.0, 0.0]

# Solve the system of equations
sol_balance = nlsolve(balance_eqs!, x0)

# Extract the harmonic balance amplitudes A and B
A_hb, B_hb = sol_balance.zero
println("Harmonic Balance Solution: A = $A_hb, B = $B_hb")

# Plot the harmonic balance solution in the time domain
tvals = 0:0.01:Tend
uhb = [A_hb * cos(omega * t) + B_hb * sin(omega * t) for t in tvals]
v_hb = [-A_hb * omega * sin(omega * t) + B_hb * omega * cos(omega * t) for t in tvals]

p5 = plot(tvals, uhb, line=:blue, title="Harmonic Balance Position", xlabel="Time", ylabel="Position")
p6 = plot(tvals, v_hb, line=:green, title="Harmonic Balance Velocity", xlabel="Time", ylabel="Velocity")
display(plot(p5, p6, layout=(1, 2)))  # Explicitly display harmonic balance results

# Compare FFT of harmonic balance with numerical solution
uf_hb = fft(uhb[1:Nsamples])
p7 = bar(fvec, 2.0 / Nsamples * log.(abs.(uf_hb[1:length(fvec)])), xlabel="Frequency", ylabel="HB Magnitude (log)", label="HB Norm")
p8 = plot(fvec, angle.(uf_hb[1:length(fvec)]), xlabel="Frequency", ylabel="HB Phase", label="HB Phase")
display(plot(p7, p8, layout=(2, 1)))  # Explicitly display FFT comparison





