using HarmonicBalance
@variables α ω ω0 F η t x(t) # declare constant variables and a function x(t)
diff_eq = DifferentialEquation(d(x,t,2) + ω0^2*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t), x)
add_harmonic!(diff_eq, x, [ω, 3ω]) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)

# solve for a range of ω
result = get_steady_states(harmonic_eq, (ω => range(0.9, 1.2, 100),
                                          α => 1., ω0 => 1.0, F => 0.01, η => 0.1))

p1=plot(result, "sqrt(u1^2 + v1^2)", legend=false)
p2=plot(result, "sqrt(u2^2 + v2^2)")
plot(p1, p2)