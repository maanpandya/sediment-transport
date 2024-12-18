using Symbolics

#t -> independent variables
#ω -> angular frequency of the ansatz
#harmonics -> list of integers for the harmonics to consider

function ansatz_definer(t, ω, harmonics)

    @variables c[1:2*length(harmonics)] #Define unknown coefficients for each harmonic term

    ansatz = zero(t)
    for i in 1:length(harmonics)
        ansatz += c[2*i-1]*sin(harmonics[i]*ω*t) + c[2*i]*cos(harmonics[i]*ω*t)
    end
    return ansatz, c 
end

#Example usage

@variables t, ω
ansatz, c = ansatz_definer(t, ω, [1, 3])
println(ansatz)
println(c)

