using SymbolicUtils:
    SymbolicUtils,
    Postwalk,
    Sym,
    BasicSymbolic,
    isterm,
    ispow,
    isadd,
    isdiv,
    ismul,
    add_with_div,
    frac_maketerm,
    @compactified,
    issym

using Symbolics:
    Symbolics,
    Num,
    unwrap,
    wrap,
    get_variables,
    Equation,
    Differential,
    @variables,
    arguments,
    simplify_fractions,
    substitute,
    term,
    expand,
    operation

using OrderedCollections: OrderedDict
using LinearAlgebra: LinearAlgebra
using SymbolicUtils
using HomotopyContinuation
using LinearAlgebra
using DynamicPolynomials
import DynamicPolynomials: coefficient

expand_all(x::Num) = Num(expand_all(x.val))
_apply_termwise(f, x::Num) = wrap(_apply_termwise(f, unwrap(x)))

"Expands using SymbolicUtils.expand and expand_exp_power (changes exp(x)^n to exp(x*n)"
function expand_all(x)
    result = Postwalk(expand_exp_power)(SymbolicUtils.expand(x))
    return isnothing(result) ? x : result
end
expand_all(x::Complex{Num}) = expand_all(x.re) + im * expand_all(x.im)

function expand_fraction(x::BasicSymbolic)
    @compactified x::BasicSymbolic begin
        Add => _apply_termwise(expand_fraction, x)
        Mul => _apply_termwise(expand_fraction, x)
        Div => sum([arg / x.den for arg in arguments(x.num)])
        _   => x
    end
end
expand_fraction(x::Num) = Num(expand_fraction(x.val))

"Apply a function f on every member of a sum or a product"
function _apply_termwise(f, x::BasicSymbolic)
    @compactified x::BasicSymbolic begin
        Add => sum([f(arg) for arg in arguments(x)])
        Mul => prod([f(arg) for arg in arguments(x)])
        Div => _apply_termwise(f, x.num) / _apply_termwise(f, x.den)
        _   => f(x)
    end
end

simplify_complex(x::Complex) = isequal(x.im, 0) ? x.re : x.re + im * x.im
simplify_complex(x) = x
function simplify_complex(x::BasicSymbolic)
    @compactified x::BasicSymbolic begin
        Add => _apply_termwise(simplify_complex, x)
        Mul => _apply_termwise(simplify_complex, x)
        Div => _apply_termwise(simplify_complex, x)
        _   => x
    end
end

"""
Perform substitutions in `rules` on `x`.
`include_derivatives=true` also includes all derivatives of the variables of the keys of `rules`.
"""
Subtype = Union{Num,Equation,BasicSymbolic}
function substitute_all(x::Subtype, rules::Dict; include_derivatives=true)
    if include_derivatives
        rules = merge(
            rules,
            Dict([Differential(var) => Differential(rules[var]) for var in keys(rules)]),
        )
    end
    return substitute(x, rules)
end
"Variable substitution - dictionary"
function substitute_all(dict::Dict, rules::Dict)::Dict
    new_keys = substitute_all.(keys(dict), rules)
    new_values = substitute_all.(values(dict), rules)
    return Dict(zip(new_keys, new_values))
end
Collections = Union{Dict,Pair,Vector,OrderedDict}
substitute_all(v::AbstractArray, rules) = [substitute_all(x, rules) for x in v]
substitute_all(x::Subtype, rules::Collections) = substitute_all(x, Dict(rules))
function substitute_all(x::Complex{Num}, rules::Collections)
    return substitute_all(x.re, rules) + im * substitute_all(x.im, rules)
end

get_independent(x::Num, t::Num) = get_independent(x.val, t)
function get_independent(x::Complex{Num}, t::Num)
    return get_independent(x.re, t) + im * get_independent(x.im, t)
end
get_independent(v::Vector{Num}, t::Num) = [get_independent(el, t) for el in v]
get_independent(x, t::Num) = x

function get_independent(x::BasicSymbolic, t::Num)
    @compactified x::BasicSymbolic begin
        Add  => sum([get_independent(arg, t) for arg in arguments(x)])
        Mul  => prod([get_independent(arg, t) for arg in arguments(x)])
        Div  => !is_function(x.den, t) ? get_independent(x.num, t) / x.den : 0
        Pow  => !is_function(x.base, t) && !is_function(x.exp, t) ? x : 0
        Term => !is_function(x, t) ? x : 0
        Sym  => !is_function(x, t) ? x : 0
        _    => x
    end
end

"Return all the terms contained in `x`"
get_all_terms(x::Num) = unique(_get_all_terms(Symbolics.expand(x).val))
function get_all_terms(x::Equation)
    return unique(cat(get_all_terms(Num(x.lhs)), get_all_terms(Num(x.rhs)); dims=1))
end
function _get_all_terms(x::BasicSymbolic)
    @compactified x::BasicSymbolic begin
        Add => vcat([_get_all_terms(term) for term in SymbolicUtils.arguments(x)]...)
        Mul => Num.(SymbolicUtils.arguments(x))
        Div => Num.([_get_all_terms(x.num)..., _get_all_terms(x.den)...])
        _   => Num(x)
    end
end
_get_all_terms(x) = Num(x)

function is_harmonic(x::Num, t::Num)::Bool
    all_terms = get_all_terms(x)
    t_terms = setdiff(all_terms, get_independent(all_terms, t))
    isempty(t_terms) && return true
    trigs = is_trig.(t_terms)

    if !prod(trigs)
        return false
    else
        powers = [max_power(first(term.val.arguments), t) for term in t_terms[trigs]]
        return all(isone, powers)
    end
end

is_harmonic(x::Equation, t::Num) = is_harmonic(x.lhs, t) && is_harmonic(x.rhs, t)
is_harmonic(x, t) = is_harmonic(Num(x), Num(t))

"Return true if `f` is a function of `var`."
is_function(f, var) = any(isequal.(get_variables(f), var))

"""
Counts the number of derivatives of a symbolic variable.
"""
function count_derivatives(x::Symbolics.BasicSymbolic)
    (Symbolics.isterm(x) || Symbolics.issym(x)) ||
        error("The input is not a single term or symbol")
    bool = Symbolics.is_derivative(x)
    return bool ? 1 + count_derivatives(first(arguments(x))) : 0
end
count_derivatives(x::Num) = count_derivatives(Symbolics.unwrap(x))


function drop_powers(expr::Num, vars::Vector{Num}, deg::Int)
    Symbolics.@variables ϵ
    subs_expr = deepcopy(expr)
    rules = Dict([var => ϵ * var for var in unique(vars)])
    subs_expr = Symbolics.expand(substitute_all(subs_expr, rules))
    max_deg = max_power(subs_expr, ϵ)
    removal = Dict([ϵ^d => Num(0) for d in deg:max_deg])
    res = substitute_all(substitute_all(subs_expr, removal), Dict(ϵ => Num(1)))
    return Symbolics.expand(res)
end

function drop_powers(expr::Vector{Num}, var::Vector{Num}, deg::Int)
    return [drop_powers(x, var, deg) for x in expr]
end

# calls the above for various types of the first argument
function drop_powers(eq::Equation, var::Vector{Num}, deg::Int)
    return drop_powers(eq.lhs, var, deg) .~ drop_powers(eq.lhs, var, deg)
end
function drop_powers(eqs::Vector{Equation}, var::Vector{Num}, deg::Int)
    return [
        Equation(drop_powers(eq.lhs, var, deg), drop_powers(eq.rhs, var, deg)) for eq in eqs
    ]
end
drop_powers(expr, var::Num, deg::Int) = drop_powers(expr, [var], deg)
drop_powers(x, vars, deg::Int) = drop_powers(Num(x), vars, deg)

"Return the highest power of `y` occurring in the term `x`."
function max_power(x::Num, y::Num)
    terms = get_all_terms(x)
    powers = power_of.(terms, y)
    return maximum(powers)
end

max_power(x::Vector{Num}, y::Num) = maximum(max_power.(x, y))
max_power(x::Complex, y::Num) = maximum(max_power.([x.re, x.im], y))
max_power(x, t) = max_power(Num(x), Num(t))

"Return the power of `y` in the term `x`"
function power_of(x::Num, y::Num)
    issym(y.val) ? nothing : error("power of " * string(y) * " is ambiguous")
    return power_of(x.val, y.val)
end

function power_of(x::BasicSymbolic, y::BasicSymbolic)
    if ispow(x) && issym(y)
        return isequal(x.base, y) ? x.exp : 0
    elseif issym(x) && issym(y)
        return isequal(x, y) ? 1 : 0
    else
        return 0
    end
end

power_of(x, y) = 0

expand_exp_power(expr::Num) = expand_exp_power(expr.val)
simplify_exp_products(x::Num) = simplify_exp_products(x.val)

"Returns true if expr is an exponential"
isexp(expr) = isterm(expr) && expr.f == exp

"Expand powers of exponential such that exp(x)^n => exp(x*n) "
function expand_exp_power(expr::BasicSymbolic)
    @compactified expr::BasicSymbolic begin
        Add => sum([expand_exp_power(arg) for arg in arguments(expr)])
        Mul => prod([expand_exp_power(arg) for arg in arguments(expr)])
        _   => ispow(expr) && isexp(expr.base) ? exp(expr.base.arguments[1] * expr.exp) : expr
    end
end
expand_exp_power(expr) = expr

"Simplify products of exponentials such that exp(a)*exp(b) => exp(a+b)
This is included in SymbolicUtils as of 17.0 but the method here avoid other simplify calls"
function simplify_exp_products(expr::BasicSymbolic)
    @compactified expr::BasicSymbolic begin
        Add => _apply_termwise(simplify_exp_products, expr)
        Div => _apply_termwise(simplify_exp_products, expr)
        Mul => simplify_exp_products_mul(expr)
        _   => expr
    end
end
function simplify_exp_products(x::Complex{Num})
    return Complex{Num}(simplify_exp_products(x.re.val), simplify_exp_products(x.im.val))
end
function simplify_exp_products_mul(expr)
    ind = findall(x -> isexp(x), arguments(expr))
    rest_ind = setdiff(1:length(arguments(expr)), ind)
    rest = isempty(rest_ind) ? 1 : prod(arguments(expr)[rest_ind])
    total = isempty(ind) ? 0 : sum(getindex.(arguments.(arguments(expr)[ind]), 1))
    if SymbolicUtils.is_literal_number(total)
        (total == 0 && return rest)
    else
        return rest * exp(total)
    end
end
simplify_exp_products(x) = x


"Expand all sin/cos powers in `x`."
function trig_reduce(x)
    x = add_div(x) # a/b + c/d = (ad + bc)/bd
    x = expand(x) # open all brackets
    x = trig_to_exp(x)
    x = expand_all(x) # expand products of exponentials
    x = simplify_exp_products(x) # simplify products of exps
    x = exp_to_trig(x)
    x = Num(simplify_complex(expand(x)))
    return simplify_fractions(x) # (a*c^2 + b*c)/c^2 = (a*c + b)/c
end

"Return true if `f` is a sin or cos."
function is_trig(f::Num)
    f = ispow(f.val) ? f.val.base : f.val
    isterm(f) && SymbolicUtils.operation(f) ∈ [cos, sin] && return true
    return false
end

"""

Returns the coefficient of cos(ωt) in `x`.
"""
function fourier_cos_term(x, ω, t)
    return _fourier_term(x, ω, t, cos)
end

"Simplify fraction a/b + c/d = (ad + bc)/bd"
add_div(x) = wrap(Postwalk(add_with_div; maketerm=frac_maketerm)(unwrap(x)))

"""

Returns the coefficient of sin(ωt) in `x`.
"""
function fourier_sin_term(x, ω, t)
    return _fourier_term(x, ω, t, sin)
end

function _fourier_term(x::Equation, ω, t, f)
    return Equation(_fourier_term(x.lhs, ω, t, f), _fourier_term(x.rhs, ω, t, f))
end

"Return the coefficient of f(ωt) in `x` where `f` is a cos or sin."
function _fourier_term(x, ω, t, f)
    term = x * f(ω * t)
    term = trig_reduce(term)
    indep = get_independent(term, t)
    ft = Num(simplify_complex(Symbolics.expand(indep)))
    ft = !isequal(ω, 0) ? 2 * ft : ft # extra factor in case ω = 0 !
    return Symbolics.expand(ft)
end

"Convert all sin/cos terms in `x` into exponentials."
function trig_to_exp(x::Num)
    all_terms = get_all_terms(x)
    trigs = filter(z -> is_trig(z), all_terms)

    rules = []
    for trig in trigs
        is_pow = ispow(trig.val) # trig is either a trig or a power of a trig
        power = is_pow ? trig.val.exp : 1
        arg = is_pow ? arguments(trig.val.base)[1] : arguments(trig.val)[1]
        type = is_pow ? operation(trig.val.base) : operation(trig.val)

        if type == cos
            term = Complex{Num}((exp(im * arg) + exp(-im * arg))^power * (1//2)^power, 0)
        elseif type == sin
            term =
                (1 * im^power) *
                Complex{Num}(((exp(-im * arg) - exp(im * arg)))^power * (1//2)^power, 0)
        end
        # avoid Complex{Num} where possible as this causes bugs
        # instead, the Nums store SymbolicUtils Complex types
        term = Num(Symbolics.expand(term.re.val + im * term.im.val))
        append!(rules, [trig => term])
    end

    result = Symbolics.substitute(x, Dict(rules))
    return convert_to_Num(result)
end
convert_to_Num(x::Complex{Num})::Num = Num(first(x.re.val.arguments))
convert_to_Num(x::Num)::Num = x

function exp_to_trig(x::BasicSymbolic)
    if isadd(x) || isdiv(x) || ismul(x)
        return _apply_termwise(exp_to_trig, x)
    elseif isterm(x) && x.f == exp
        arg = first(x.arguments)
        trigarg = Symbolics.expand(-im * arg) # the argument of the to-be trig function
        trigarg = simplify_complex(trigarg)

        # put arguments of trigs into a standard form such that sin(x) = -sin(-x), cos(x) = cos(-x) are recognized
        if isadd(trigarg)
            first_symbol = minimum(
                cat(string.(arguments(trigarg)), string.(arguments(-trigarg)); dims=1)
            )

            # put trigarg => -trigarg the lowest alphabetic argument of trigarg is lower than that of -trigarg
            # this is a meaningless key but gives unique signs to all sums
            is_first = minimum(string.(arguments(trigarg))) == first_symbol
            return if is_first
                cos(-trigarg) - im * sin(-trigarg)
            else
                cos(trigarg) + im * sin(trigarg)
            end
        end
        return if ismul(trigarg) && trigarg.coeff < 0
            cos(-trigarg) - im * sin(-trigarg)
        else
            cos(trigarg) + im * sin(trigarg)
        end
    else
        return x
    end
end

exp_to_trig(x) = x
exp_to_trig(x::Num) = exp_to_trig(x.val)
exp_to_trig(x::Complex{Num}) = exp_to_trig(x.re) + im * exp_to_trig(x.im)


function ansatz_simplifier(ansatz)
    return trig_reduce(ansatz)
end






#t -> independent variables
#ω -> angular frequency of the ansatz
#harmonics -> list of integers for the harmonics to consider
#c -> list of unknown coefficients

function ansatz_definer(t, ω, harmonics, n)

    @variables c[1:2*length(harmonics)*n] #Define unknown coefficients for each harmonic term

    ansatz = [] #Each element is the ansatz for one variable
    global k = 1
    for j in 1:n
        push!(ansatz, zero(t))
        for i in 1:length(harmonics)
            ansatz[j] += c[(2*k-1)]*sin(harmonics[i]*ω*t) + c[(2*k)]*cos(harmonics[i]*ω*t)
            k += 1
        end
    end
    
    return ansatz, c 
end

#ansatz -> truncated trigonometric series from ansatz_definer function
#powers -> list of the exponents which will be used to calculate the powers of the ansatz
#derivatives -> list of the orders of the derivatives to calculate for the ansatz

function power_derivatives(t, ansatz, powers, derivatives)

    ansatz_derivatives = []
    ansatz_powers = []
    for k in 1:length(ansatz)
        #Derivative calculation
        sub_ansatz_derivatives = []
        i = 1
        while i <= length(derivatives)
            
            if i == 1
                push!(sub_ansatz_derivatives, Symbolics.simplify(Symbolics.derivative(ansatz[k], t))) 
            else 
                push!(sub_ansatz_derivatives, Symbolics.simplify(Symbolics.derivative(sub_ansatz_derivatives[i-1], t)))
            end
            i += 1
        end

        #Power calculation

        sub_ansatz_powers = []
        for j in 1:length(powers)
            push!(sub_ansatz_powers, ansatz_simplifier(ansatz[k]^(powers[j])))
        end
        push!(ansatz_derivatives, sub_ansatz_derivatives)
        push!(ansatz_powers, sub_ansatz_powers)
    end

    return ansatz_powers, ansatz_derivatives
end


"Return the coefficient of f(ωt) in `x` where `f` is a cos or sin."
function _fourier_term(x, ω, t, f)
    term = x * f(ω * t)
    term = trig_reduce(term)
    indep = get_independent(term, t)
    ft = Num(simplify_complex(Symbolics.expand(indep)))
    ft = !isequal(ω, 0) ? 2 * ft : ft # extra factor in case ω = 0 !
    return Symbolics.expand(ft)
end

# Harmonic Balance Substitution
function harmonic_balance_substitution(ansatz, ode, ansatz_powers, ansatz_derivatives, harmonics, truncation_level)
    # Simplification Rules for Trigonometric Identities
    r1 = @rule cos((~x))^2 => 0.5 + 0.5*cos(2*(~x))
    r2 = @rule sin((~x))^2 => 0.5 - 0.5*sin(2*(~x))
    r3 = @rule 2.0*~a*~b*sin((~x))*cos((~x)) => ~a*~b*sin(2*(~x))
    r4 = @rule 2.0*~a*~b*cos((~x))*sin((~x)) => ~a*~b*sin(2*(~x))
    r5 = @rule cos((~x))^3 => 0.75*cos((~x)) + 0.25*cos(3*(~x)) 
    r6 = @rule sin((~x))^3 => 0.75*sin((~x)) - 0.25*sin(3*(~x)) 
    ruleset = RuleSet([r1, r2, r3, r4, r5, r6])

    # Combine ansatz, powers, and derivatives into a substitution dictionary
    combined_dict = Dict(
        ansatz => ansatz,
        (ansatz)^2 => ansatz_powers[1],
        (ansatz)^3 => ansatz_powers[2],
        Differential(t)(ansatz) => ansatz_derivatives[1],
        Differential(t)(Differential(t)(ansatz)) => ansatz_derivatives[2]
    )

    # Step 1: Explicit substitution
    substituted_eq = Symbolics.substitute(ode, combined_dict)
    println(" After Substitution: ", substituted_eq)

    # Step 2: Expand and simplify the equation
    expanded_eq = Symbolics.expand(substituted_eq)
    println(" After Expansion: ", expanded_eq)
    
    simplified_eq = Symbolics.simplify(expanded_eq, ruleset)
    println(" After Simplification: ", simplified_eq)

    # Step 3: Define valid harmonics (up to truncation level)
    valid_harmonics = [sin(n * ω * t) for n in harmonics] ∪ [cos(n * ω * t) for n in harmonics]
    #println("\n Valid Harmonics: ", valid_harmonics)

    # Step 4: Group terms by harmonics
    #grouped_terms = Dict(harmonic => Num(0) for harmonic in valid_harmonics)
    terms = isa(simplified_eq, Symbolics.Add) ? Symbolics.arguments(simplified_eq) : [simplified_eq]
    grouped_terms = Dict(harmonic => Num(0) for harmonic in valid_harmonics)
    #println(" Terms after Expansion: ", terms)


    for term in terms
        term_expr = term isa Equation ? term.lhs : term
        term_expr = Symbolics.expand(term_expr)

        matched = false
        for harmonic in valid_harmonics
            try
                # Match terms explicitly to harmonics
                if occursin(string(harmonic), string(term_expr))
                    grouped_terms[harmonic] += term_expr
                    matched = true
                    break
                end
            catch e
                @warn " Failed to match term: $term_expr with harmonic: $harmonic. Error: $e"
            end
        end

        if !matched
            @warn " Unmatched term during harmonic grouping: $term_expr"
        end
    end

    # Step 5: Display grouped terms
    #println(" Grouped Terms: ", grouped_terms)

    # Step 6: Generate harmonic balance equations
    harmonic_equations = [grouped_terms[harmonic] ~ 0 for harmonic in valid_harmonics if !iszero(grouped_terms[harmonic])]

    #println(" Harmonic Equations: ", harmonic_equations)
    return harmonic_equations
end


function harmonic_separation_with_fourier(equations::Vector{Equation}, ω, t)
    harmonics = [
        (sin, ω, t),
        (cos, ω, t)#,
        # (sin, 3*ω, t),
        # (cos, 3*ω, t)
    ]

    # Initialize an empty list to store harmonic coefficients
    harmonic_coefficients = []

    for eq in equations
        lhs = eq.lhs  # Extract left-hand side of the current equation

        for (f, freq, time) in harmonics
            coeff = _fourier_term(lhs, freq, time, f)
            harmonic = f(freq * time)
            push!(harmonic_coefficients, (harmonic, coeff))  # Store as tuple (harmonic, coefficient)
        end
    end

    return harmonic_coefficients
end

#Example usage
harmonics = [1]
NumMasses = 2
@variables t ω δ α β γ c[1:2*length(harmonics)*NumMasses]
# Use this once we make the functiosn array friendly
# if NumMasses == 1
#     @variables t ω δ α β γ c[1:2*length(harmonics)*NumMasses]
# else
#     @variables t ω δ[1:NumMasses] α[1:NumMasses] β[1:NumMasses] γ m[1:NumMasses] c[1:2*length(harmonics)*NumMasses]
# end
D = Differential(t)
ansatz, c = ansatz_definer(t, ω, harmonics, NumMasses)
println(ansatz)

#duffing_eq = D(D(ansatz)) + δ*D(ansatz) + α*ansatz + β*(ansatz)^3 ~ γ*cos(ω*t)
duffing_eq = [
    D(D(ansatz[1])) ~ -δ*D(ansatz[1]) - α*ansatz[1] - β*(ansatz[1])^3 + δ*(D(ansatz[2])-D(ansatz[1])) + α*(ansatz[2]-ansatz[1]) + β*(ansatz[2]-ansatz[1])^3,
    D(D(ansatz[2])) ~ γ*cos(ω*t) - δ*(D(ansatz[2])-D(ansatz[1])) - α*(ansatz[2]-ansatz[1]) - β*(ansatz[2] - ansatz[1])^3
]
println(duffing_eq)
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

println("The harmonic balance substitution is:")
harmonic_equations = harmonic_balance_substitution(ansatz, duffing_eq, ansatz_powers, ansatz_derivatives, [1], 1)
println("Output of harmonic balance substitution:")
println(harmonic_equations)

# Perform harmonic separation and store coefficients in a list
harmonic_coefficients = harmonic_separation_with_fourier(harmonic_equations, ω, t)
# Display harmonic coefficients
for (harmonic, coeff) in harmonic_coefficients
    println("Harmonic: $harmonic, Coefficient: $coeff")
end
println(harmonic_coefficients)

input_funcs = [coeff for (harmonic, coeff) in harmonic_coefficients]
println(input_funcs)


function solve_polynomial_system(num_harmonics,input_alpha, input_beta, input_gamma, input_delta, input_omega, input_funcs, returnnonsingular=false)
    n = num_harmonics * 2 # number of unknown variables

    # 1. Declare polynomial variables
    @polyvar c_p[1:n]  # coefficients
    @polyvar α_poly[1:length(input_alpha)]
    @polyvar β_poly[1:length(input_beta)]
    @polyvar γ_poly[1:length(input_gamma)]
    @polyvar δ_poly[1:length(input_delta)]
    @polyvar ω_poly[1:length(input_omega)]
    @variables γ[1:length(input_gamma)]

    # add gamma as forcing term to first equation from input_funcs
    input_funcs[2] -= γ[1]

    # 2. Create substitution dictionary
    substitution_dict = Dict()
    
    # Map symbolic variables to polynomial variables
    for i in 1:length(input_alpha)
        substitution_dict[α[i]] = input_alpha[i]
        substitution_dict[β[i]] = input_beta[i]
        substitution_dict[γ[i]] = input_gamma[i]
        substitution_dict[δ[i]] = input_delta[i]
    end
    
    for i in 1:length(input_omega)
        substitution_dict[ω[i]] = input_omega[i]
    end

    # Map coefficients
    for i in 1:n
        substitution_dict[c[i]] = c_p[i]  # Map symbolic c to polynomial c
    end

    # 3. Convert each equation in input_funcs to polynomial form
    poly_funcs = [Symbolics.substitute(eq, substitution_dict) for eq in input_funcs]

    println("Converted polynomial system:")
    println(poly_funcs)
    println(typeof(poly_funcs))
    
    # Define the system of equations
    f = poly_funcs
    println("The system of equations is:")
    println(f)

    # Define the system as a polynomial system
    system = HomotopyContinuation.System(f)

    # Solve the system
    result = HomotopyContinuation.solve(system)

    # Collect real solutions
    real_solutions = []
    for sol in result
        if HomotopyContinuation.is_real(sol)
            push!(real_solutions, sol)
        end
    end

    # Return nonsingular solutions if requested
    if returnnonsingular
        nonsingular_solutions = []
        for sol in real_solutions
            if HomotopyContinuation.is_nonsingular(sol)
                push!(nonsingular_solutions, sol)
            end
        end
        return real_solutions, nonsingular_solutions
    end

    return real_solutions
end

println("The results are:")
println(solve_polynomial_system(2,[1], [0.04], [1], [0.1], [1], input_funcs))