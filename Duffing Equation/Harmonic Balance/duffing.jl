using HarmonicBalance
@variables α, ω, ω0, F, t, γ, x(t); # declare constant variables and a function x(t)
# define ODE 
diff_eq = DifferentialEquation(d(x,t,2) + ω0^2*x + α*x^3 + γ*d(x,t) ~ F*cos(ω*t), x)

# specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)
add_harmonic!(diff_eq, x, ω)

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = (α => 1., ω0 => 1.0, F => 0.01, γ=>0.01) # fixed parameters
varied = ω => LinRange(0.9, 1.2, 100)              # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)

plot(result, "sqrt(u1^2 + v1^2)", size=(400,250))