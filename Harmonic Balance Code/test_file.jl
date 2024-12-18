using Symbolics
using SymbolicUtils

@variables k, A, B

r1 = @rule cos(~x)^3 => 0.75 * cos(~x) + 0.25 * cos(3 * ~x)
r2 = @rule sin(~x)^3 => 0.75 * sin(~x) - 0.25 * sin(3 * ~x)
r3 = @rule cos(~x)^2 => 1 - sin(~x)^2
r4 = @rule sin(~x)^2 => 1 - cos(~x)^2
r5 = @rule sin(3*~x) => 0 
r6 = @rule cos(3*~x) => 0 

expr1 = simplify(expand((A*cos(k)+B*sin(k))^2), RuleSet([r1,r2,r3,r4]))
expr2 = simplify(expand(expr1), RuleSet([r1,r2,r3,r4]))
expr3 = simplify(expand(expr2), RuleSet([r5,r6]))
