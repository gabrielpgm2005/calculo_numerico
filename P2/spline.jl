#VICTOR DE OLIVEIRA SILVA E GABRFIEL PERROUT GOMES DE MOURA

function spline_natural(xk :: Vector{<:Real}, yk :: Vector{<:Real})
    n = length(xk)
    @assert n == length(yk) "Os vetores xk e yk devem ter o mesmo tamanho."
    @assert n >= 2 "São necessários pelo menos 2 pontos para construir o spline."
    @assert issorted(xk) "Os pontos xk devem estar em ordem crescente."

    T = float(promote_type(eltype(xk), eltype(yk)))
    m = n - 1  # número de subintervalos
    h = [T(xk[i+1] - xk[i]) for i in 1:m]

    # Sistema linear A * b_full = rhs para b_1, ..., b_n
    A = zeros(T, n, n)
    rhs = zeros(T, n)

    # Condições de contorno naturais: S''(x_1) = 0 e S''(x_n) = 0
    A[1, 1] = one(T)
    A[n, n] = one(T)

    # Continuidade de S' e S'' nos nós internos
    for i in 1:(n - 2)
        A[i+1, i]   = h[i]
        A[i+1, i+1] = 2 * (h[i] + h[i+1])
        A[i+1, i+2] = h[i+1]
        rhs[i+1]    = 3 * ((yk[i+2] - yk[i+1]) / h[i+1] - (yk[i+1] - yk[i]) / h[i])
    end

    b_full = A \ rhs

    a = Vector{T}(undef, m)
    b = b_full[1:m]
    c = Vector{T}(undef, m)
    d = Vector{T}(undef, m)

    for i in 1:m
        a[i] = (b_full[i+1] - b_full[i]) / (3 * h[i])
        c[i] = (yk[i+1] - yk[i]) / h[i] - (h[i] / 3) * (b_full[i+1] + 2 * b_full[i])
        d[i] = T(yk[i])
    end

    return a, b, c, d
end

function avalia_spline(
    xk :: Vector{<:Real},
    a :: Vector{<:Real},
    b :: Vector{<:Real},
    c :: Vector{<:Real},
    d :: Vector{<:Real},
    x :: Real
)
    m = length(a)
    @assert length(xk) == m + 1 && length(b) == m && length(c) == m && length(d) == m "Dimensões incompatíveis entre xk e os coeficientes."

    # Busca binária pelo subintervalo [xk[i], xk[i+1]] que contém x
    i = clamp(searchsortedlast(xk, x), 1, m)
    dx = x - xk[i]

    # Avaliação pelo método de Horner
    return ((a[i] * dx + b[i]) * dx + c[i]) * dx + d[i]
end

xk = [0, 1, 3]
yk = [1, 2, 0]

a, b, c, d = spline_natural(xk, yk)
println("a = ", a) # [-0.3333333333333333, 0.16666666666666666]
println("b = ", b) # [0.0, -1.0]
println("c = ", c) # [1.3333333333333333, 0.33333333333333326]
println("d = ", d) # [1.0, 2.0]

println("S(0.5) = ", avalia_spline(xk, a, b, c, d, 0.5)) # 1.625
println("S(2)   = ", avalia_spline(xk, a, b, c, d, 2))   # 1.5
