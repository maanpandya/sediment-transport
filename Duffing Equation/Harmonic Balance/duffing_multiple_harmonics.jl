using HarmonicBalance

n_cases = 100 #Defines the number of different driving frequencies to try
harmonics = [1] #Defines the harmonic to use as an integer which multiplies the fundamental harmonic
print_harmonic_eqs = true #Set to print or not the harmonic equations

#α -> linear stiffness coefficient
#γ -> amplitude of periodic driving force
#β -> non-linear stiffness coefficient
#δ -> damping coefficient
#ω -> angular frequency of periodic driving force  

@variables α, ω, β, δ, t, γ, x(t); # declare constant variables and a function x(t) 

diff_eq = DifferentialEquation(d(x,t,2) + α*x + β*x^3 + δ*d(x,t) ~ γ*cos(ω*t), x)

fixed = (α => 1.0, β => 0.04, δ => 0.1, γ=>1)   #fixed parameters
varied = ω => range(0.01, 3, n_cases)           #range of driving frequencies

for i in 1:length(harmonics)
    add_harmonic!(diff_eq, x, [harmonics[i]*ω]) #add the ith-harmonic
end

harmonic_eq = get_harmonic_equations(diff_eq)

if print_harmonic_eqs
    println(harmonic_eq)
end

result = get_steady_states(harmonic_eq, varied, fixed)
println(result)

#Branch 1 always contains the real and stable solutions
stable_classification = classify_branch(result, "stable")
physical_classification = classify_branch(result, "physical")
println("Stable and physical solution are located in the following branches:")
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

    for g in 1:2:2*length(harmonics)
        local u = sol[key_list[g]]
        local v = sol[key_list[g+1]]
        push!(single_sol_coeffs, u)
        push!(single_sol_coeffs, v)
    end
    local ω_case = real(sol[key_list[2*length(harmonics)+1]])
    push!(single_sol_coeffs, ω_case)
    push!(sol_coeffs, single_sol_coeffs)
end

println(sol_coeffs)



#plot(result, "sqrt(u1^2 + v1^2)")
#savefig("C:\\Users\\LENOVO\\Desktop\\duffing_multiple_response.png")

#if plot_harmonic_coefficients
#    plot(xlabel="ω", ylabel="u" title="Harmonic Coefficients", size=(700, 600), titlepadding=100)
#
#end