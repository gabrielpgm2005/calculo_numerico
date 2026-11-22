#VICTOR DE OLIVEIRA SILVA E GABRFIEL PERROUT GOMES DE MOURA

function lagrange(xk :: Vector{<:Real}, yk :: Vector{<:Real}, x :: Real)
    n = length(xk)
    @assert n == length(yk) "Os vetores xk e yk devem ter o mesmo tamanho."

    T = float(promote_type(eltype(xk), eltype(yk), typeof(x)))
    P = zero(T)

    for i in 1:n
        L = one(T)
        for j in 1:n
            if j != i
                L *= (x - xk[j]) / (xk[i] - xk[j])
            end
        end
        P += yk[i] * L
    end

    return P
end

xk = [1, 4, 9]
yk = [1, 2, 3]
x = 2

println(lagrange(xk, yk, x)) # 1.3666666666666667