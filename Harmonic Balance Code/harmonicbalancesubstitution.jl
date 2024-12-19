using Symbolics
using SymbolicUtils

# t -> independent variables
# ω -> angular frequency of the ansatz
# harmonics -> list of integers for the harmonics to consider
# c -> list of unknown coefficients

function ansatz_definer(t, ω, harmonics)
    @variables c[1:2*length(harmonics)] #Define unknown coefficients for each harmonic term

    ansatz = zero(t)
    for i in 1:length(harmonics)
        ansatz += c[2*i-1]*sin(harmonics[i]*ω*t) + c[2*i]*cos(harmonics[i]*ω*t)
    end
    return ansatz, c 
end

# ansatz -> truncated trigonometric series from ansatz_definer function
# powers -> list of the exponents which will be used to calculate the powers of the ansatz
# derivatives -> list of the orders of the derivatives to calculate for the ansatz

function power_derivatives(t, ansatz, powers, derivatives)
    # Derivative calculation
    ansatz_derivatives = []
    i = 1
    while i <= length(derivatives)
        if i == 1
            push!(ansatz_derivatives, Symbolics.simplify(Symbolics.derivative(ansatz, t))) 
        else 
            push!(ansatz_derivatives, Symbolics.simplify(Symbolics.derivative(ansatz_derivatives[i-1], t)))
        end
        i += 1
    end

    # Power calculation
    r1 = @rule cos((~x))^2 => 0.5 + 0.5*cos(2*(~x))
    r2 = @rule sin((~x))^2 => 0.5 - 0.5*sin(2*(~x))
    r3 = @rule 2.0*~a*~b*sin((~x))cos((~x)) => ~a~b*sin(2*(~x))
    r4 = @rule 2.0*~a*~b*cos((~x))sin((~x)) => ~a~b*sin(2*(~x))
    r5 = @rule cos((~x))^3 => 0.75*cos((~x)) + 0.25*cos(3*(~x)) 
    r6 = @rule sin((~x))^3 => 0.75*sin((~x)) - 0.25*sin(3*(~x)) 

    ansatz_powers = []
    for j in 1:length(powers)
        power_expr = ansatz^(powers[j])
        power_simpl = simplify(expand(power_expr), RuleSet([r1, r2, r3, r4, r5, r6]))
        push!(ansatz_powers, power_simpl)
    end

    return ansatz_powers, ansatz_derivatives
end

function harmonic_balance_substitution(ansatz, ode, powers_derivatives, harmonics, truncation_level)
    """
    Substitutes the ansatz into the ODE, removes higher harmonics, and groups terms.

    Parameters:
        ansatz: Truncated Fourier series ansatz.
        ode: Symbolic differential equation.
        powers_derivatives: Dictionary mapping powers and derivatives for substitution.
        harmonics: List of base harmonics (e.g., [sin(ωt), cos(ωt), sin(2ωt), cos(2ωt)]).
        truncation_level: Truncation level (max harmonic n).

    Returns:
        List of polynomial equations grouped by harmonics.
    """

    # Rules for simplification to remove higher harmonics
    r1 = @rule cos((~x))^2 => 0.5 + 0.5*cos(2*(~x))
    r2 = @rule sin((~x))^2 => 0.5 - 0.5*sin(2*(~x))
    r3 = @rule 2.0*~a*~b*sin((~x))cos((~x)) => ~a~b*sin(2*(~x))
    r4 = @rule 2.0*~a*~b*cos((~x))sin((~x)) => ~a~b*sin(2*(~x))
    r5 = @rule cos((~x))^3 => 0.75*cos((~x)) + 0.25*cos(3*(~x)) 
    r6 = @rule sin((~x))^3 => 0.75*sin((~x)) - 0.25*sin(3*(~x)) 

    ruleset = RuleSet([r1, r2, r3, r4, r5, r6])

    # Substitute powers and derivatives into the ODE
    substituted_eq = substitute(ode, powers_derivatives)
    expanded_eq = expand(substituted_eq)  # expand the substituted equation
    simplified_eq = simplify(expand(substituted_eq), ruleset)  # Simplify using the rules

    # Filter harmonics to keep only terms up to truncation_level
    filtered_harmonics = [sin(n*ω*t) for n in 1:truncation_level] ∪ [cos(n*ω*t) for n in 1:truncation_level]

    # Group terms by harmonics and extract coefficients
    grouped_terms = Dict()
    for harmonic in filtered_harmonics
        #grouped_terms[harmonic] = coefficient(expanded_eq, harmonic)  # extract coefficient of harmonic
        grouped_terms[harmonic] = coefficient(simplified_eq, harmonic)
    
    # Form polynomial equations by equating coefficients of each harmonic to zero
    polynomial_eqs = [grouped_terms[harmonic] ~ 0 for harmonic in filtered_harmonics]
    
    return polynomial_eqs
end

# Example usage
@variables t ω δ α β γ F u₁ v₁
D = Differential(t)

# Define parameters
δ = 0.1    # damping coefficient
ω = 1.5    # frequency  
α = 1.0    # linear stiffness
β = 0.04   # nonlinearity
γ = 1.0    # driving force amplitude
n = 2      # number of unknown variables
harmonics = [1] 
truncation_level = 1

# Define the ODE(eg)
duffing_eq = D(D(u₁*cos(ω*t) + v₁*sin(ω*t))) + δ*D(u₁*cos(ω*t) + v₁*sin(ω*t)) + α*(u₁*cos(ω*t) + v₁*sin(ω*t)) + β*(u₁*cos(ω*t) + v₁*sin(ω*t))^3 ~ F*cos(ω*t)

ansatz, c = ansatz_definer(t, ω, harmonics)
println(ansatz)

ansatz_powers, ansatz_derivatives = power_derivatives(t, ansatz, [2, 3], [1, 2])

# Define harmonics to consider for harmonic balance
harmonics_to_consider = [sin(ω * t), cos(ω * t)]

polynomial_eqs = harmonic_balance_substitution(ansatz, duffing_eq, ansatz_derivatives, harmonics_to_consider, truncation_level)


println("Polynomial equations:")
for eq in polynomial_eqs
    println(eq)
    end
end