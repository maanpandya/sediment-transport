using HarmonicBalance

n_cases = 100 #Defines the number of different driving frequencies to try
n_harmonics = 1 #Defines the number of harmonics to use
print_harmonic_eqs = true #Set to print or not the harmonic equations

#α -> linear stiffness coefficient
#γ -> amplitude of periodic driving force
#β -> non-linear stiffness coefficient
#δ -> damping coefficient
#ω -> angular frequency of periodic driving force  

@variables α, ω, β, δ, t, γ, x(t); # declare constant variables and a function x(t) 

diff_eq = DifferentialEquation(d(x,t,2) + α*x + β*x^3 + δ*d(x,t) ~ γ*cos(ω*t), x)

fixed = (α => 1.0, β => 1.0, δ => 0.3, γ=>0.3)   #fixed parameters
varied = ω => range(1.0, 1.5, n_cases)           #range of driving frequencies

harmonics = filter(x -> x % 2 != 0, collect(1:2*n_harmonics)) #Harmonics or multiples of ω
for i in harmonics
    add_harmonic!(diff_eq, x, [i*ω]) #add the ith-harmonic
end

harmonic_eq = get_harmonic_equations(diff_eq)

if print_harmonic_eqs
    println(harmonic_eq)
end

result = get_steady_states(harmonic_eq, varied, fixed)
print(result)
#Branch 1 always contains the real and stable solutions
stable_classification = classify_branch(result, "stable")
physical_classification = classify_branch(result, "physical")
println(stable_classification)
println(physical_classification)

#sol = get_single_solution(result, branch=1, index=2)
#println(sol)

sol_coeffs = [] #List will contain the coefficients for all n_cases
for j in 1:n_cases
    local sol = get_single_solution(result, branch=1, index=j) #Solution for 1 one of the driving frequency values
    local single_sol_coeffs = [] #List to contain the Fourier coefficients for this particular solution
    #Access the contents and keys of the solution dictionary variable
    local key_list = []
    for (key, value) in sol
        #println("Key: ", key, ", Value: ", value)
        push!(key_list, key)
    end

    for g in 1:n_harmonics
        local u = real(sol[key_list[harmonics[g]]])
        local v = real(sol[key_list[harmonics[g]+1]])
        push!(single_sol_coeffs, u)
        push!(single_sol_coeffs, v)
    end
    push!(sol_coeffs, single_sol_coeffs)
end
println(sol_coeffs)
plot(result, "sqrt(u1^2 + v1^2)", title="1st Harmonic (Mass 1)", size=(700, 600), legend=false)
