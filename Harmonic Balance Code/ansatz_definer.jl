using Symbolics

#t -> independent variables
#ω -> angular frequency of the ansatz
#harmonics -> list of integers for the harmonics to consider

function ansatz_definer(t, ω, harmonics)

    ansatz = 0
    for i in 1:length(harmonics)
        ansatz = ansatz + sin(harmonics[i]*ω*t) + cos(harmonics[i]*ω*t)
    end
    return ansatz 
end

#Example usage
"""
@variables t, ω
ansatz = ansatz_definer(t, ω, [1, 3])
println(ansatz)
"""