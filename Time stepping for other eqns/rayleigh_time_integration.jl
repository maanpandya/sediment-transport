using DifferentialEquations
using Plots

#Define the Duffing equations


function rayleigh(du, u, p, t)
    ε, γ, ω = p
    du[1] = u[2]
    du[2] = -u[1] - ε * (u[2]^3 - u[2]) + γ * cos(ω * t)
end

# Setting parameters
ε = 0.8
γ = 2.0
ω = 1.0
p = [ε, γ, ω]

# Initial conditions, time span, solver steps stay the same
u0 = [0.0, 0.0]
tspan = (0.0, 150.0)
dt = 0.001
# Solve the ODE
prob = ODEProblem(rayleigh, u0, tspan, p)
sol1 = solve(prob, ImplicitEuler(), dt=dt)
sol2 = solve(prob, Euler(), dt=dt)

# Plot for both solvers and compare
plot(sol1, vars=(0, 1), label="Implicit Scheme")
plot!(sol2, vars=(0, 1), label="Explicit Scheme")

# import libary for fft
using FFTW
using Statistics

# performing fft on the solution u1 to see the frequency components

Nsamples = 2000
Tstart = 100
Tend = 150
dt = (Tend - Tstart) / (Nsamples)
tsampled = Vector(Tstart:dt:Tend)
usampled = [sol1(t)[1] for t in tsampled]

# performing fft
fsampled = fft(usampled)
# fsampled = fsampled / Nsamples


# frequency
f_max = 1 / (2 * dt)
fstep = 1/((Tend - Tstart))
fvec = Vector(0:fstep:f_max)

# determining the dominant frequency
fmax = fvec[argmax(abs.(fsampled[1:length(fvec)]))]
println("Dominant harmonic: ", fmax*2*π)

# find amplitude of the dominant frequency
# Find index of maximum amplitude in FFT data
max_idx = argmax(abs.(fsampled[1:length(fvec)]))

# Calculate amplitude (with correct scaling)
dominant_amplitude = 2.0/Nsamples * abs(fsampled[max_idx])

# println("Dominant frequency: ", fvec[max_idx], " rad/s")
println("Amplitude at dominant frequency: ", dominant_amplitude)


#keep horizonal axis as ω
fvec = [2 * π * f for f in fvec]
# onlt plot the first half of the frequency components
fvec = fvec[1:Int(Nsamples/20)]
fsampled = fsampled[1:Int(Nsamples/20)]

# only plot first 50 harmonics
fvec = fvec[1:60]
fsampled = fsampled[1:60]


# plotting the frequency components
# p1 = plot(fvec, 2.0/Nsamples * abs.(fsampled[1:length(fvec)]), label="FFT of u1", xlabel="ω", ylabel="Amplitude")
p2 = plot(fvec, angle.(fsampled[1:length(fvec)]), label="FFT of u1", xlabel="Frequency", ylabel="Phase")
# plot(p1, p2, layout=(2,1))

p1 = plot(fvec, 2.0/Nsamples * abs.(fsampled[1:length(fvec)]),
    label="FFT of u1",
    xlabel="ω",
    ylabel="Amplitude",
    xticks=0:1:maximum(fvec),  # Set ticks every 1 unit
    # Alternative: for even denser ticks use
    # xticks=0:0.5:maximum(fvec),  # Set ticks every 0.5 units
    grid=false
)


# show the plot
plot(p1)

# save the plot
savefig("rayleigh_fft_time.png")