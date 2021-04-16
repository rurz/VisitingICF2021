# Dependences
using HypergeometricFunctions
using SpecialFunctions

export kp, binom, Φ

"`kp(n, q, N)` gives the symmetric Kravchuk polynomiar for n ∈[0, N] and q ∈ [-N,N]"
kp(n, q, N) = HypergeometricFunctions._₂F₁general2(BigFloat(-n), BigFloat(-q), -N, 2.0)

A(n, j) = ((-1.0)^n) / (2.0^j)

binom(n, k) = (1/n) * (1/beta(k + 1, n - k + 1))

B(n, q, j) = √(binom(2 * j, n) * binom(2 * j, j + q))

"Fundamental definition of the Kravchuk function"
ϕ(n::Integer, q::Integer, j::Integer) = A(n, j) * B(n, q, j) * kp(n, j + q, 2 * j)

"`Φ(n, q, j)` gives the Kravchuk function with the ranges n ∈ [0, 2j], q ∈ [-j, j]. It is tested to give correct values for j ≤ 128."
function Φ(n, q, j)
    if 0 ≤ n ≤ j
        return Float64(ϕ(n, q, j))
    elseif j < n ≤ 2 * j
        return Float64(((-1.0)^q) * ϕ(2 * j - n, q, j))
    end
end
