# Used dependences
using WignerFunctions

"This script evaluates the Laguerre-Kravchuk functions"

"`tp(l, m₁, m₂)` is the redefinition of the tupled of indexes for the Wigner little-_d_ matrix"
tp(l, m₁, m₂) = WignerFunctions.Index((l, m₁, m₂))

"`ld(β, l, m₁, m₂) is the redefinition of the Wigner little-_d_ matrix`"
ld(β, l, m₁, m₂) = WignerFunctions.naive(β, tp(l, m₁, m₂))

"`λ(n, j₁, j₂) is a romboid counter for the little-_d_ matrix`"
function λ(n, j₁, j₂)
    if 0 ≤ n ≤ 2 * j₂
        return (1/2) * n
    elseif 2 * j₂ < n ≤ 2 * j₁
        return j₂
    elseif 2 * j₁ ≤ n ≤ 2 * (j₁ + j₂)
        return j₁ + j₂ - (1/2) * n
    end
end

"`μ(n, j₁, j₂, n₁, n₂) is a romboid counter for the little-_d_ matrix`"
function μ(n, j₁, j₂, n₁, n₂)
    if 0 ≤ n ≤ 2 * j₂
        return (1/2) * (n₁ - n₂)
    elseif 2 * j₂ < n ≤ 2 * j₁
        return j₂ - n₂
    elseif 2 * j₁ ≤ n ≤ 2 * (j₁ + j₂)
        return (1/2) * (n₁ - n₂) - j₁ + j₂
    end
end

expf(σ, n₁, n₂) = exp(σ * 1im * (1/4) * π * (n₁ - n₂))
