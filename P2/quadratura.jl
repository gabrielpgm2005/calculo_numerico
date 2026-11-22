#VICTOR DE OLIVEIRA SILVA E GABRFIEL PERROUT GOMES DE MOURA

function pontos_pesos_gauss(::Type{T}, n :: Int) where {T <: AbstractFloat}
    t = zeros(T, n)
    w = zeros(T, n)
    m = fld(n + 1, 2) # número de raízes não-negativas

    for i in 1:m
        # Aproximação inicial de Tricomi para a i-ésima raiz positiva
        ti = cos(T(pi) * (4 * i - 1) / (4 * n + 2))

        # Refinamento pelo Método de Newton usando a recorrência de Legendre
        for _ in 1:50
            p0, p1 = one(T), ti
            dp1 = one(T)

            for k in 2:n
                pk = ((2 * k - 1) * ti * p1 - (k - 1) * p0) / k
                dpk = k * p1 + ti * dp1
                p0, p1 = p1, pk
                dp1 = dpk
            end

            delta = p1 / dp1
            ti -= delta
            if abs(delta) <= eps(T) * max(one(T), abs(ti))
                break
            end
        end

        # Avalia P'_n(ti) no ponto convergido para obter o peso w_i
        p0, p1 = one(T), ti
        dp1 = one(T)
        for k in 2:n
            pk = ((2 * k - 1) * ti * p1 - (k - 1) * p0) / k
            dpk = k * p1 + ti * dp1
            p0, p1 = p1, pk
            dp1 = dpk
        end

        wi = 2 / ((1 - ti^2) * dp1^2)

        # Preenche simetricamente em [-1, 1]
        t[n - i + 1] = ti
        w[n - i + 1] = wi
        t[i] = -ti
        w[i] = wi
    end

    if isodd(n)
        t[m] = zero(T)
    end

    return t, w
end

function gauss(f :: Function, a :: Real, b :: Real, n :: Int)
    @assert 2 <= n <= 10 "O número de pontos n deve estar entre 2 e 10."

    T = float(promote_type(typeof(a), typeof(b)))
    t, w = pontos_pesos_gauss(T, n)

    meio = (T(a) + T(b)) / 2
    metade_compr = (T(b) - T(a)) / 2

    soma = zero(T)
    for i in 1:n
        x = meio + metade_compr * t[i]
        soma += w[i] * f(x)
    end

    return metade_compr * soma
end

# Teste 1: f(x) = e^x em [0, 1] (Exato: e - 1 ≈ 1.7182818285)
println("Teste 1 (e^x):")
println("n = 2: ", gauss(exp, 0, 1, 2))   # 1.7178963780314205
println("n = 3: ", gauss(exp, 0, 1, 3))   # 1.718281004372522
println("n = 5: ", gauss(exp, 0, 1, 5))   # 1.718281828458938

# Teste 2: f(x) = √x em [0, 1] (Exato: 2/3 ≈ 0.6666666667)
println("\nTeste 2 (√x):")
println("n = 2:  ", gauss(sqrt, 0, 1, 2))  # 0.6738873386790492
println("n = 5:  ", gauss(sqrt, 0, 1, 5))  # 0.6672967897192363
println("n = 10: ", gauss(sqrt, 0, 1, 10)) # 0.6667560428791446