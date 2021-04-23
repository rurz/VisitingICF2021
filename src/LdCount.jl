"Construction of the little-d Wigner function, and the counting modes for rotations"

using WignerFunctions

export tp, ld, λ, μ

"`tp(l, m₁, m₂)` is the redefinition of the tupled of indexes for the Wigner little-_d_ matrix"
tp(l, m₁, m₂) = WignerFunctions.Index((l, m₁, m₂))

"`ld(β, l, m₁, m₂) is the redefinition of the Wigner little-_d_ matrix`"
ld(β, l, m₁, m₂) = WignerFunctions.naive(β, tp(l, m₁, m₂))

"`λ(n, j₁, j₂)` is a romboid counter for the little-_d_ matrix"
function λ(n, j₁, j₂)
    if 0 ≤ n ≤ 2 * j₂
        return (1/2) * n
    elseif 2 * j₂ < n ≤ 2 * j₁
        return j₂
    elseif 2 * j₁ ≤ n ≤ 2 * (j₁ + j₂)
        return j₁ + j₂ - (1/2) * n
    end
end

"`μ(n, m, j₁, j₂)` is a romboid counter for the little-_d_ matrix. We have the equivalence m = n₁ - n₂, and n = n₁ + n₂"
function μ(n, m, j₁, j₂)
    if 0 ≤ n ≤ 2 * j₂
        return (1/2) * m
    elseif 2 * j₂ < n ≤ 2 * j₁
        return j₂ + (1/2) * (m - n)
    elseif 2 * j₁ ≤ n ≤ 2 * (j₁ + j₂)
        return (1/2) * m - j₁ + j₂
    end
end

"`expf(σ, m)` is a index-dependent phase"
expf(σ, n₁, n₂) = exp(σ * 1im * (1/4) * π * m)

function Ld(j₁, j₂)
    v₀ = zeros((2 * j₁ + 1, 2 * j₂ + 1))
    for n in 0:2 * (j₁ + j₂)
        for m in -n:2:n
            v₀[n + 1, m + n + 1] = ld(π/2, λ(n, j₁, j₂), μ(n, m, j₁, j₂), μ(n, m, j₁, j₂))
        end
    end
    return v₀
end
