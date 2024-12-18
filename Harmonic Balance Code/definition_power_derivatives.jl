using Symbolics
using SymbolicUtils

#t -> independent variables
#ω -> angular frequency of the ansatz
#harmonics -> list of integers for the harmonics to consider
#c -> list of unknown coefficients

function ansatz_definer(t, ω, harmonics)

    @variables c[1:2*length(harmonics)] #Define unknown coefficients for each harmonic term

    ansatz = zero(t)
    for i in 1:length(harmonics)
        ansatz += c[2*i-1]*sin(harmonics[i]*ω*t) + c[2*i]*cos(harmonics[i]*ω*t)
    end
    return ansatz, c 
end

#ansatz -> truncated trigonometric series from ansatz_definer function
#powers -> list of the exponents which will be used to calculate the powers of the ansatz
#derivatives -> list of the orders of the derivatives to calculate for the ansatz

function power_derivatives(t, ansatz, powers, derivatives)

    #Derivative calculation
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

    #Power calculation

    r1 = @rule cos((~x))^2 => 0.5 + 0.5*cos(2*(~x))
    r2 = @rule sin((~x))^2 => 0.5 - 0.5*sin(2*(~x))
    r3 = @rule 2.0*~a*~b*sin((~x))*cos((~x)) => ~a*~b*sin(2*(~x))
    r4 = @rule 2.0*~a*~b*cos((~x))*sin((~x)) => ~a*~b*sin(2*(~x))
    r5 = @rule cos((~x))^3 => 0.75*cos((~x)) + 0.25*cos(3*(~x)) 
    r6 = @rule sin((~x))^3 => 0.75*sin((~x)) - 0.25*sin(3*(~x)) 

    ansatz_powers = []
    for j in 1:length(powers)
        local power = simplify(expand(ansatz^(powers[j])), RuleSet([r1, r2, r3, r4, r5, r6]))
        push!(ansatz_powers, simplify(expand(power), RuleSet([r1, r2, r3, r4, r5, r6])))
    end

    return ansatz_powers, ansatz_derivatives
end

#Example usage
@variables t, ω
ansatz, c = ansatz_definer(t, ω, [1])
println(ansatz)

ansatz_powers, ansatz_derivatives = power_derivatives(t, ansatz, [2, 3], [1, 2])
println("The results are:")
println("Power of 2")
println(ansatz_powers[1])
println("Power of 3")
println(ansatz_powers[2])
println("First derivative")
println(ansatz_derivatives[1])
println("Second derivative")
println(ansatz_derivatives[2])
