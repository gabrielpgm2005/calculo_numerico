#VICTOR DE OLIVEIRA SILVA E GABRFIEL PERROUT GOMES DE MOURA

function ponto_medio(f :: Function, a :: Real, b :: Real, m :: Int)
    @assert m >= 1 "O número de subintervalos m deve ser positivo."

    T = float(promote_type(typeof(a), typeof(b)))
    h = (T(b) - T(a)) / m
    soma = zero(T)

    for i in 1:m
        x_mid = T(a) + (i - T(0.5)) * h
        soma += f(x_mid)
    end

    return h * soma
end

function trapezios(y :: Vector{<:Real}, h :: Real)
    n = length(y)
    @assert n >= 2 "O vetor y deve conter pelo menos 2 pontos."

    T = float(promote_type(eltype(y), typeof(h)))
    soma = T(y[1]) + T(y[end])

    for i in 2:(n - 1)
        soma += 2 * y[i]
    end

    return (T(h) / 2) * soma
end

function trapezios(f :: Function, a :: Real, b :: Real, m :: Int)
    @assert m >= 1 "O número de subintervalos m deve ser positivo."

    T = float(promote_type(typeof(a), typeof(b)))
    h = (T(b) - T(a)) / m
    y = [f(T(a) + i * h) for i in 0:m]

    return trapezios(y, h)
end

function simpson(y :: Vector{<:Real}, h :: Real)
    m = length(y) - 1
    @assert m >= 2 && iseven(m) "O número de subintervalos (length(y) - 1) deve ser par e >= 2."

    T = float(promote_type(eltype(y), typeof(h)))
    soma = T(y[1]) + T(y[end])

    for i in 2:m
        peso = iseven(i) ? 4 : 2
        soma += peso * y[i]
    end

    return (T(h) / 3) * soma
end

function simpson(f :: Function, a :: Real, b :: Real, m :: Int)
    @assert m >= 2 && iseven(m) "O número de subintervalos m deve ser par e >= 2."

    T = float(promote_type(typeof(a), typeof(b)))
    h = (T(b) - T(a)) / m
    y = [f(T(a) + i * h) for i in 0:m]

    return simpson(y, h)
end

f(x) = 1 / x

println("Ponto médio: ", ponto_medio(f, 1, 2, 4)) # 0.6912198912198912
println("Trapézios:   ", trapezios(f, 1, 2, 4))   # 0.6970238095238095
println("Simpson:     ", simpson(f, 1, 2, 4))     # 0.6932539682539683
println("Exato (ln 2):", log(2))                  # 0.6931471805599453

# Teste com pontos tabelados
y = [1, 1/1.25, 1/1.5, 1/1.75, 1/2]
h = 0.25

println("Trapézios (y): ", trapezios(y, h))       # 0.6970238095238095
println("Simpson (y):   ", simpson(y, h))         # 0.6932539682539683