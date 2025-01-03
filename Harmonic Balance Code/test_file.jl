using Symbolics
using SymbolicUtils

@syms x::Real A::Real B::Real C::Real w::Real t::Real y::Real
 
r1 = @acrule sin(~w * ~t)cos(~w * ~t) => 0.5*sin(2* ~w * ~t)
r2 = @acrule cos(~w * ~t)^2 => 0.5 + 0.5*cos(2* ~w * ~t)
r3 = @acrule sin(~w * ~t)^2 => 0.5 - 0.5*cos(2* ~w * ~t)
r4 = @acrule sin(2* ~w * ~t) => 0 
r5 = @acrule cos(2* ~w * ~t) => 0
r6 = @acrule cos(~w * ~t)sin(~w * ~t) => ~a* ~b*0.5*sin(2* ~w * ~t)
r7 = @acrule cos(~w * ~t)cos(2* ~w * ~t) => 0.5*cos(~w * ~t) + 0.5*cos(3* ~w * ~t)
r8 = @acrule sin(~x)cos(~y) => 0.5*(sin(~x + ~y) + sin(~x - ~y))
r12 = @acrule sin((~-)*(~a * ~w * ~t)) => - sin((~a * ~w * ~t))

expr20 = simplify(expand(A*sin(2*w*t)cos(13*w*t)), RuleSet([r8, r12]))
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

