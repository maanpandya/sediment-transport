using Symbolics
using SymbolicUtils

@syms x::Real A::Real B::Real C::Real w::Real t::Real y::Real
 
r1 = @acrule sin(~w * ~t)cos(~w * ~t) => 0.5*sin(2* ~w * ~t)
r2 = @acrule cos(~w * ~t)^2 => 0.5 + 0.5*cos(2* ~w * ~t)
r3 = @acrule sin(~w * ~t)^2 => 0.5 - 0.5*cos(2* ~w * ~t)
r4 = @acrule sin(2* ~w * ~t) => 0 
r5 = @acrule cos(2* ~w * ~t) => 0
r8 = @acrule sin(~x)cos(~y) => 0.5*(sin(~x + ~y) + sin(~x - ~y))
r12 = @acrule sin((~-)*(~a * ~w * ~t)) => -1*sin((~a * ~w * ~t))

expr20 = simplify(expand(A*sin(2*w*t)cos(13*w*t)), RuleSet([r8]))
println(expr20)
expr21 = simplify(expand(expr20), RuleSet([r12]))
println(expr21)
"""
expr1 = simplify(expand((A*cos(w*t)+B*sin(w*t))^2))
println(expr1)
expr2 = simplify(expand(expr1), RuleSet([r1]))
println(expr2)
expr3 = simplify(expand(expr2), RuleSet([r2,r3]))
println(expr3)
expr4 = simplify(expand(expr3), RuleSet([r4,r5]))
println(expr4)
"""

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

