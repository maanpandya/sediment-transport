using HarmonicBalance
using Symbolics
@variables α, ω, ω0, F, t, γ, x(t); # declare constant variables and a function x(t)
# define ODE 
diff_eq = DifferentialEquation(d(x,t,2) + ω0^2*x + α*x^3 + γ*d(x,t) ~ F*cos(ω*t), x)

fixed = (α => 10., ω0 => 3, F => 5, γ=>0.01)   # fixed parameters
varied = ω => range(0.9, 1.4, 2)           # range of parameter values

add_harmonic!(diff_eq, x, [ω]) # specify the two-harmonics ansatz
harmonic_eq = get_harmonic_equations(diff_eq)
print(harmonic_eq)

result = get_steady_states(harmonic_eq, varied, fixed)

print(result)
#print(size(result))
#println(get_variables(result))
sol = get_single_solution(result, branch=2, index=2)
#print(sol)
#print(keys(sol))  # Replace `solution` with your OrderedDict variable

key_list = []
for (key, value) in sol
    println("Key: ", key, ", Value: ", value)
    push!(key_list, key)
end
println(key_list)
u1_value = sol[key_list[1]]
v1_value = sol[key_list[2]]

println("The ansatz coefficients are:")
println("u1: $u1_value")
println("v1: $v1_value")



#print(result[1])

#print(typeof(get_single_solution(result; branch=1, index=(1))))

#p1=plot(result, "sqrt(u1^2 + v1^2)", legend=false, filename =" coupled_duffing_sols ")
#p2=plot(result, "sqrt(u2^2 + v2^2)")
#plot(p1, p2, size=(800,350))