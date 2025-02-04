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
    issym,
    simplify

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

using Plots
using OrderedCollections: OrderedDict
using LinearAlgebra: LinearAlgebra
using SymbolicUtils
using HomotopyContinuation
using LinearAlgebra, SparseArrays
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
    # Simplification Rules
    r1 = @rule cos((~x))^2 => 0.5 + 0.5*cos(2*(~x))
    r2 = @rule sin((~x))^2 => 0.5 - 0.5*sin(2*(~x))
    r3 = @rule 2.0*~a*~b*sin((~x))*cos((~x)) => ~a*~b*sin(2*(~x))
    r4 = @rule 2.0*~a*~b*cos((~x))*sin((~x)) => ~a*~b*sin(2*(~x))
    r5 = @rule cos((~x))^3 => 0.75*cos((~x)) + 0.25*cos(3*(~x)) 
    r6 = @rule sin((~x))^3 => 0.75*sin((~x)) - 0.25*sin(3*(~x)) 
    ruleset = RuleSet([r1, r2, r3, r4, r5, r6])

    # Ensure ODE is a vector of equations
    if !(ode isa Vector)
        ode = [ode]
    end

    substituted_eqs = []
    for eq in ode #This way, each ODE sees all the needed substitutions (including the other masses) and you won’t be left with stray Differential(t)(...).
        # Build substitution dictionary for ansatz i
        combined_dict = Dict()
        for (idx, ans) in pairs(ansatz)
            combined_dict[ans] = ans
            combined_dict[ans^2] = ansatz_powers[idx][1]
            combined_dict[ans^3] = ansatz_powers[idx][2]
            combined_dict[Differential(t)(ans)] = ansatz_derivatives[idx][1]
            combined_dict[Differential(t)(Differential(t)(ans))] = ansatz_derivatives[idx][2]
        end
        # (Optionally merge in other ansatz[j] substitutions if needed)

        # Substitution
        substituted_eq = Symbolics.substitute(eq, combined_dict)
        expanded_eq = Symbolics.expand(substituted_eq)
        simplified_eq = Symbolics.simplify(expanded_eq, ruleset)

        push!(substituted_eqs, simplified_eq)
    end

    # Define valid harmonics
    valid_harmonics = [sin(n * ω * t) for n in harmonics] ∪ [cos(n * ω * t) for n in harmonics]

    all_harmonic_equations = []
    for simplified_eq in substituted_eqs
        # Combine LHS and RHS => eq_eval = LHS - RHS
        eq_eval = simplified_eq.lhs - simplified_eq.rhs
        eq_eval = Symbolics.expand(eq_eval)
        eq_eval = Symbolics.simplify(eq_eval, ruleset)

        # Collect terms
        terms = isa(eq_eval, Symbolics.Add) ? Symbolics.arguments(eq_eval) : [eq_eval]
        grouped_terms = Dict(harmonic => Num(0) for harmonic in valid_harmonics)

        for term in terms
            matched = false
            for harmonic in valid_harmonics
                if occursin(string(harmonic), string(term))
                    grouped_terms[harmonic] += term
                    matched = true
                    break
                end
            end
            if !matched
                @warn "Unmatched term: $term"
            end
        end

        # Generate equations
        for harmonic in valid_harmonics
            if !iszero(grouped_terms[harmonic])
                push!(all_harmonic_equations, grouped_terms[harmonic] ~ 0)
            end
        end
    end

    return all_harmonic_equations
end


function harmonic_separation_with_fourier(equations, ω, t, length_harmonics)
    if length_harmonics == 1
        harmonics = [
            (sin, ω, t),
            (cos, ω, t)
        ]
    elseif length_harmonics == 2
        harmonics = [
            (sin, ω, t),
            (cos, ω, t),
            (sin, 3*ω, t),
            (cos, 3*ω, t)
        ]
    end

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

function create_diff_matrix(nx, L=2π; disc, bc)
    dx = L / (nx - 1)
    A = zeros(nx, nx)
    if disc == :central
        for i in 2:nx-1
            A[i, i-1] = -1/(2dx)
            A[i, i+1] = 1/(2dx)
        end
    elseif disc == :upwind
        for i in 2:nx-1
            A[i, i-1] = -1/dx
            A[i, i] = 1/dx
        end
    elseif disc == :downwind
        for i in 2:nx-1
            A[i, i] = -1/dx
            A[i, i+1] = 1/dx
        end
    end
    # Periodic boundary conditions
    if bc == :periodic
        A[1, end] = -1/(2dx)
        A[1, 2] = 1/(2dx)
        A[end, end-1] = -1/(2dx)
        A[end, 1] = 1/(2dx)
    elseif bc == :nonperiodic
        A[1, 1] = -1/dx
        A[1, 2] = 1/dx
        A[end, end-1] = -1/dx
        A[end, end] = 1/dx
    elseif bc == :upwind
        A[1, end] = -1/dx
        A[1, 1] = 1/dx
        A[end, end-1] = -1/dx
        A[end, end] = 1/dx
    elseif bc == :downwind
        A[1, 1] = -1/dx
        A[1, 2] = 1/dx
        A[end, end] = -1/dx
        A[end, 1] = 1/dx
    end
    return A
end

"""
The root issue is that even when using central (or downwind) differences 
the advection‐equation (a hyperbolic PDE) has little numerical diffusion 
so that high–frequency (or aliasing) errors appear as spurious 
discontinuities (here near π/2). One common fix is to add a small amount 
of artificial viscosity that damps those oscillations.
"""

function create_dd2_matrix(nx, L=2π; disc=:central, bc=:periodic)
    dx = L/(nx-1)
    D2 = zeros(nx, nx)
    if disc == :central
        for i in 2:nx-1
            D2[i,i-1] = 1/(dx^2)
            D2[i,i]   = -2/(dx^2)
            D2[i,i+1] = 1/(dx^2)
        end
    end
    if bc == :periodic
        D2[1,end]   = 1/(dx^2)
        D2[1,1]     = -2/(dx^2)
        D2[1,2]     = 1/(dx^2)
        D2[end,end-1] = 1/(dx^2)
        D2[end,end]   = -2/(dx^2)
        D2[end,1]     = 1/(dx^2)
    end
    return D2
end

#Example usage
harmonics = [1, 3]
nx = 1000
ωtest = 1.0
disc = :central
bc = :periodic
A = create_diff_matrix(nx, disc=disc ,bc=bc)
#print(A)
dx = 2π/(nx-1)
ν = dx/10
B = create_dd2_matrix(nx, 2π; disc=:central, bc=:periodic)
@variables t ω F[1:nx]
# F is just ones for now
F = ones(nx)
# Explicitly define forcing terms without broadcasting
f = [F[i] * sin(ω*t) for i in 1:nx]
n = nx
D = Differential(t)
ansatz, c = ansatz_definer(t, ω, harmonics, n)
println("Size of c " * string(size(c)))
if nx < 4
    println("The ansatz is:")
    println(ansatz)
end

# Time derivatives of the ansatz
du = [Symbolics.derivative(u,t) for u in ansatz]
if nx < 4
    println(du)
end

ansatz_powers, ansatz_derivatives = power_derivatives(t, ansatz, [2, 3], [1, 2])
println("ansatz powers and derivatives done")

"""
This modification adds a small dissipative term that damps out the high–frequency aliasing responsible for the discontinuity at some points. 
Adjust ν until the spurious jump is removed without overly damping the true advection.
"""

# Discretized PDE: du/dt = A*u + f
ode_system = Equation[]
for i in 1:nx
    lhs = du[i]
    rhs = sum(A[i,j]*ansatz[j] for j in 1:nx) + ν*sum(B[i,j]*ansatz[j] for j in 1:nx) + f[i]
    push!(ode_system, lhs ~ rhs)
end
#println(ode_system)

println("The harmonic balance substitution is:")
harmonic_equations = harmonic_balance_substitution(ansatz, ode_system, ansatz_powers, ansatz_derivatives, harmonics, 1)
println("Output of harmonic balance substitution:")
#println(harmonic_equations)

# Perform harmonic separation and store coefficients in a list
harmonic_coefficients = harmonic_separation_with_fourier(harmonic_equations, ω, t, length(harmonics))
# Display harmonic coefficients
# for (harmonic, coeff) in harmonic_coefficients
#     println("Harmonic: $harmonic, Coefficient: $coeff")
# end
#println(harmonic_coefficients)

input_funcs = [coeff for (harmonic, coeff) in harmonic_coefficients]
println("The input functions are:")
#println(input_funcs)


function solve_polynomial_system(n, input_alpha, input_beta, input_gamma, input_delta, input_omega, input_funcs, returnnonsingular=false)

    # 1. Declare polynomial variables
    @polyvar c_p[1:n]  # coefficients

    # 2. Create substitution dictionary
    substitution_dict = Dict()
    
    # Map symbolic variables to polynomial variables
    for (i, val) in enumerate(input_alpha)
        substitution_dict[α] = val
    end
    for (i, val) in enumerate(input_beta)
        substitution_dict[β] = val
    end
    for (i, val) in enumerate(input_gamma)
        substitution_dict[γ] = val
    end
    for (i, val) in enumerate(input_delta)
        substitution_dict[δ] = val
    end
    
    for (i, val) in enumerate(input_omega)
        substitution_dict[ω] = val
    end

    # Map coefficients
    for i in eachindex(c_p)
        substitution_dict[c[i]] = c_p[i]  # Map symbolic c to polynomial c
    end

    # 3. Convert each equation in input_funcs to polynomial form
    poly_funcs = [Symbolics.substitute(eq, substitution_dict) for eq in input_funcs]

    println("Converted polynomial system:")
    println(poly_funcs)
    println(length(poly_funcs))
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

# function solve_harmonic_balance(input_funcs, c, ω_value=1.0)
#     # To solve a system of polynomial equations
# end

function extract_matrix_from_equations(input_funcs, n)
    # n is the length of vector c (6 in this case)
    H = zeros(n, n)
    b = zeros(n)
    
    for (i, eq) in enumerate(input_funcs)
        # Convert equation to string and split into terms
        eq_str = string(eq)
        terms = split(eq_str, r"(?=[-+])")
        
        for term in terms
            term = strip(term)
            if term == ""
                continue
            end
            
            # Handle constant terms
            if !contains(term, "c[") && !contains(term, "ω")
                b[i] = parse(Float64, term)
                continue
            end
            
            # Extract coefficient and index for c terms
            if contains(term, "c[")
                coef = 1.0
                if term[1] == '-'
                    coef = -1.0
                end
                
                # Extract numerical coefficient if present
                m = match(r"([-+]?[\d.]+)?c\[(\d+)\]", term)
                if m !== nothing
                    if m.captures[1] !== nothing
                        coef *= parse(Float64, m.captures[1])
                    end
                    idx = parse(Int, m.captures[2])
                    H[i, idx] = coef
                end
            end
            
            # Handle ω terms
            if contains(term, "ω")
                coef = 1.0
                if term[1] == '-'
                    coef = -1.0
                end
                
                # Extract which c[i] multiplies ω
                m = match(r"c\[(\d+)\]\*ω", term)
                if m !== nothing
                    idx = parse(Int, m.captures[1])
                    H[i, idx] = coef
                end
            end
        end
    end
    
    return H, b
end

println("The results are:")
H, b = extract_matrix_from_equations(input_funcs, 2*length(harmonics)*nx)
println("Matrix H:")
#display(H)
println("\nVector b:")
#display(b)

# Solve the system
c_p = H \ (-b)

println("\nSolved coefficients:")
#display(c_p)

"""
Now because nx = 3, it means we discretized space into 3 points. The behavior of the PDE at each point we assume is given by the ansatz. If harmonics = [1],
then we are only considering the first harmonic. This means that the ansatz will have 2 terms for each point. This is why we have 6 coefficients in the system
and 6 equations, one for each coefficient. c_p is the solution to the system of equations and it gives us the coefficients for the ansatz.
So right now, the first 2 values of c_p are the coefficients for the first point (sine and cosine), the next 2 values are for the second point and the last 2 values are for the third point.
Because we have these coefficients, we can now substitute them back into the ansatz to get the solution to the PDE at each point. And then we can plot the solution with this grid.
"""

# Create substitution dictionary for coefficients c using c_p 
substitution_dict = Dict{Any,Any}()
for i in 1:length(c)
    substitution_dict[c[i]] = Num(c_p[i])
end

#Set value of ω
substitution_dict[ω] = Num(ωtest)

# Substitute the coefficients back into the ansatz for each spatial point
u_points = [Symbolics.substitute(u, substitution_dict) for u in ansatz]

# Evaluate and print u(x,t) at each spatial grid point: x = 2π*k/(nx-1), for k = 0 to nx-1
# println("Substituted ansatz solutions for u(x,t):")
# for k in 0:(nx - 1)
#     x = 2π * k/(nx - 1)
#     println("u(x=$(x), t) = ", u_points[k+1])
# end

"""
So in the case of 3 spatial points, these are the solutions to the PDE at each point. We can now plot these solutions to see how the PDE behaves at each point.
u(x=0.0, t) = -0.9480228105904878cos(t) + 0.3265823128063383sin(t)
u(x=3.141592653589793, t) = -0.9972983717810876cos(t) - 0.016974830730531663sin(t)
u(x=6.283185307179586, t) = -1.0546788176284247cos(t) + 0.3096074820758066sin(t)
"""

# Plot the solutions
# Evaluate u(x,t) at each spatial grid point for time values in t_values
x_values = [2π * k/(nx - 1) for k in 0:(nx - 1)]
step = 0.1
t_values = collect(0:step:10)

# u_values is a matrix with dimensions (nx, length(t_values))
u_values = [Symbolics.substitute(u, t => t_val) for u in u_points, t_val in t_values]
# println(size(u_values))
# println(typeof(u_values[1][1]))
# println(u_values)

# Convert the symbolic expressions to numerical values
# Adjust conversion as needed; here we use Symbolics.evaluate to obtain a Float64
u_numeric = [Float64(Symbolics.unwrap(expr)) for expr in u_values]
print("okay")
# println(u_numeric)
# println(size(u_numeric))
# The contour function expects the z matrix arranged with rows matching the t_values.
# Create the contour plot
p = contour(t_values, x_values, u_numeric, xlabel="Time", ylabel="Space", title="Contour plot of u(x,t), nx = $nx, $disc", fill=true)
umax = maximum(u_numeric)
# Save the plot
savefig(p, "C:/Python Code/sediment-transport/Scalar Wave Equation/Viscous1Dscalarcontour.png")

# Create an animated plot: space on x-axis, u on y-axis, time evolves
anim = @animate for t0 in t_values
    # For each time t0, evaluate u at each spatial point
    u_at_t = [Float64(Symbolics.unwrap(Symbolics.substitute(u, t => t0))) for u in u_points]
    plot(x_values, u_at_t,
         xlabel="Space", ylabel="u(x,t)",
         title="Time: $(round(t0, digits=2)), nx = $nx, $disc",
         legend=false, ylim=(-umax*1.25,umax*1.25))
end

gif(anim, "C:/Python Code/sediment-transport/Scalar Wave Equation/Viscous1Dscalarsolution.gif", fps=30)