using Symbolics, HarmonicBalance, OrderedCollections

N = 5 #Number of intervals
h = 1/N #Step size
n_cases = 5 #Number of driving frequencies to test

@variables t #Independent variable
@variables (u(t))[1:N+1]  #Node variables
@variables A, k, ω, c #Physical variables

#Create and store the spatially discretized PDE
discretized_ODEs = Vector{Equation}() 

for i in 1:N+1
    if i == 1 || i == N+1 #Left or right boundary
        push!(discretized_ODEs, d(u[i], t) ~ A*sin(k*u[i] - ω*t))
    
    else
        push!(discretized_ODEs, d(u[i], t) ~ A*sin(k*u[i] - ω*t) - (c/(2*h))*(u[i+1] - u[i-1]))
        
    end
end

nodes = collect(u)  #Convert Symbolics.Arr to Vector{Num}
diff_eq = DifferentialEquation(discretized_ODEs, nodes) #Create the DifferentialEquation object

#Define parameter values
fixed = (A => 2, c => 10, k => 62) #Fixed parameters
varied = ω => range(0.01, 3, n_cases) #Range of driving frequencies

[add_harmonic!(diff_eq,node, [ω]) for node in nodes]

println(diff_eq)
harmonic_eq = get_harmonic_equations(diff_eq)
println(harmonic_eq)