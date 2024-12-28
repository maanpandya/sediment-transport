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
    
    r1 = @acrule cos((~x))^2 => 0.5 + 0.5*cos(2*(~x))
    r2 = @acrule sin((~x))^2 => 0.5 - 0.5*cos(2*(~x))

    r3 = @acrule sin((~x))cos((~x)) => sin(2*(~x))*0.5
    r4 = @acrule cos((~x))sin((~x)) => sin(2*(~x))*0.5

    r5 = @acrule sin(~ω * ~t)cos(2* ~ω * ~t) => 0.5*sin(3* ~ω * ~t) - 0.5*sin(~ω * ~t)
    r6 = @acrule cos(~t * ~ω)cos(2* ~t * ~ω) => 0.5*cos(~t * ~ω) + 0.5*cos(3* ~t * ~ω) 

    r7 = @acrule cos((~x))^3 => 0.75*cos((~x)) + 0.25*cos(3*(~x)) 
    r8 = @acrule sin((~x))^3 => 0.75*sin((~x)) - 0.25*sin(3*(~x))

    r9 = @acrule sin(~x)cos(~y) => 0.5*(sin(~x + ~y) + sin(~x - ~y))
    r10 = @acrule sin(~x)sin(~y) => 0.5*(cos(~x - ~y) - cos(~x + ~y))
    r11 = @acrule cos(~x)cos(~y) => 0.5*(cos(~x - ~y) + cos(~x + ~y))

    r12 = @acrule sin(-1 * (~t * ~ω)) => -1 * sin(~t * ~ω)
    r13 = @acrule cos(-1 * (~t * ~ω)) => cos(~t * ~ω)

    ansatz_powers = []
    for j in 1:length(powers)
        local power1 = simplify(expand(ansatz^(powers[j])))
        #println(power1)
        local power9 = simplify(expand(power1), RuleSet([r1, r2]))
        #println(power9)
        local power2 = simplify(expand(power9), RuleSet([r3, r4]))
        #println(power2)
        local power3 = simplify(expand(power2), RuleSet([r7, r8])) 
        #println(power3)
        local power4 = simplify(expand(power3), RuleSet([r9, r10, r11]))
        #println(power4)
        local power5 = simplify(expand(power4), RuleSet([r12, r13]))
        push!(ansatz_powers, simplify(expand(power5)))
    end

    return ansatz_powers, ansatz_derivatives
end

#Example usage
@variables t, ω
ansatz, c = ansatz_definer(t, ω, [1, 2])
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


