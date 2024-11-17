using HarmonicBalance

n_cases = 100 #Defines the number of different driving frequencies to try
harmonics = [1, 3] #Defines the harmonic to use as an integer which multiplies the fundamental harmonic
print_harmonic_eqs = true #Set to print or not the harmonic equations

#α1 -> linear stiffness coefficient of first spring
#β1 -> non-linear stiffness coefficient of first spring
#δ1 -> damping coefficient of first spring
#α2 -> linear stiffness coefficient of second spring
#β2 -> non-linear stiffness coefficient of second spring
#δ2 -> damping coefficient of second spring
#ω -> angular frequency of periodic driving force 
#γ -> amplitude of periodic driving force 
#m1 -> mass of the first mass
#m2 -> mass of the second mass

@variables α1, β1, δ1, α2, β2, δ2, m1, m2, t, γ, x1(t), x2(t); # declare constant variables and a function x(t) 

free_eq = [m1*d(x1,t,2) + δ1*d(x1,t) + α1*x1 + β1*x1^3 - δ2*(d(x2,t) - d(x1,t)) - α2*(x2-x1) - β2*(x2-x1)^3,
m2*d(x2,t,2) + δ2*(d(x2,t) - d(x1,t)) + α2*(x2-x1) + β2*(x2-x1)^3]

forces = [0, γ*cos(ω*t)]

diff_eq = DifferentialEquation(free_eq - forces, [x1, x2])

fixed = (α1 => 1.0, β1 => 1.0, δ1 => 0.1, α2 => 1.0, β2 => 1.0, δ2 => 0.1, γ=>0.37, m1=>1, m2=>1)   #fixed parameters
varied = ω => range(0.01, 3, n_cases)           #range of driving frequencies

for i in 1:length(harmonics)
    add_harmonic!(diff_eq, x1, [harmonics[i]*ω]) #add the ith-harmonic
    add_harmonic!(diff_eq, x2, [harmonics[i]*ω]) #add the ith-harmonic
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
    println(sol)
    local single_sol_coeffs = [] #List to contain the Fourier coefficients for this particular solution
    #Access the contents and keys of the solution dictionary variable
    local key_list = []
    for (key, value) in sol
        #println("Key: ", key, ", Value: ", value)
        push!(key_list, key)
    end

    for g in 1:2:4*length(harmonics)
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
#xlim=(0, 3), ylim=(0, 9)
p1 = plot(result, "sqrt(u1^2 + v1^2)", title="1st Harmonic (Mass 1)", size=(700, 600), legend=false)
p2 = plot(result, "sqrt(u2^2 + v2^2)", title="2nd Harmonic (Mass 1)", size=(700, 600), legend=false)
p3 = plot(result, "sqrt(u3^2 + v3^2)", title="1st Harmonic (Mass 2)", size=(700, 600), legend=false)
p4 = plot(result, "sqrt(u4^2 + v4^2)", title="2nd Harmonic (Mass 2)", size=(700, 600))
plot(p1, p2, p3, p4)
savefig("C:\\Users\\LENOVO\\Desktop\\duffing_multiple_response.png")

#if plot_harmonic_coefficients
#    plot(xlabel="ω", ylabel="u" title="Harmonic Coefficients", size=(700, 600), titlepadding=100)
#
#end