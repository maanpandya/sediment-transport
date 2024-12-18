using Symbolics, LinearAlgebra, SparseArrays, Plots

# Define the constants for the PDE
const D = 1.0  # Diffusion coefficient
const N = 32    # Grid size (N x N)
const _DD = 100.0  # Scale factor for the diffusion term

# Create the spatial grid
const X = reshape([i for i in 1:N for j in 1:N], N, N)
const Y = reshape([j for i in 1:N for j in 1:N], N, N)

# Discretize the second-order spatial derivatives using finite differences
const Mx = Array(Tridiagonal([1.0 for i in 1:N-1], [-2.0 for i in 1:N], [1.0 for i in 1:N-1]))
const My = copy(Mx)

# Apply boundary conditions for the finite difference operator
Mx[2, 1] = 2.0
Mx[end-1,end] = 2.0
My[1, 2] = 2.0
My[end,end-1] = 2.0

# Define the discretized diffusion equation as an ODE function
function f(u, p, t)
    # u is a 3D array where u[:,:,1] is the state at time t
    A = u[:,:,1]
    
    # Apply the diffusion operator in both x and y directions
    MyA = My * A
    AMx = A * Mx
    DA = @. _DD * (MyA + AMx)  # Diffusion term
    
    # The rate of change (time derivative) of the concentration
    dA = @. D * DA
    
    return dA
end

@variables u[1:N, 1:N, 1]
dA = simplify.(f(collect(u), nothing, 0.0))


println(dA[3, 3])



