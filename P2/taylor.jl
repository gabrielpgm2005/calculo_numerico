#VICTOR DE OLIVEIRA SILVA E GABRFIEL PERROUT GOMES DE MOURA

function taylor(d :: Vector{<:Real}, x0 :: Real, x :: Real)
    @assert !isempty(d) "O vetor d deve conter pelo menos o valor de f(x0)."

    T = float(promote_type(eltype(d), typeof(x0), typeof(x)))
    n = length(d) - 1
    h = T(x) - T(x0)

    # Avaliação pelo esquema de Horner adaptado à série de Taylor
    soma = T(d[end])
    for k in n:-1:1
        soma = T(d[k]) + (h / k) * soma
    end

    return soma
end

# Teste 1: f(x) = e^x em torno de x0 = 0 (grau 4)
d_exp = [1, 1, 1, 1, 1]
println("e^1 (Taylor grau 4):   ", taylor(d_exp, 0, 1))   # 2.7083333333333335
println("e^1 (exato):           ", exp(1))               # 2.718281828459045

# Teste 2: f(x) = sen(x) em torno de x0 = 0 (grau 5)
d_sin = [0, 1, 0, -1, 0, 1]
println("sen(0.5) (Taylor g.5): ", taylor(d_sin, 0, 0.5)) # 0.47942708333333334
println("sen(0.5) (exato):      ", sin(0.5))             # 0.479425538604203