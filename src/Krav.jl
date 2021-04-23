"Construction of the Kravchuk functions for at most 256 points (j ≤ 128)"

using HypergeometricFunctions
using SpecialFunctions

export kp, binom, Φ

"`kp(n, q, N)` gives the symmetric Kravchuk polynomia for n ∈ [0, N] and q ∈ [-N,N]"
kp(n, q, N) = HypergeometricFunctions._₂F₁general2(-BigFloat(n), -q, -N, 2)

A(n, j) = ((-1.0)^n) / (2.0^j)

binom(n, k) = (1/n) * (1/beta(k + 1, n - k + 1))

B(n, q, j) = √(binom(2 * j, n) * binom(2 * j, j + q))

"Fundamental definition of the Kravchuk function"
ϕ(n, q, j) = Float64(A(n, j) * B(n, q, j) * kp(n, j + q, 2 * j))

"`Φ(n, q, j)` gives the Kravchuk function with the ranges n ∈ [0, 2j], q ∈ [-j, j]. It is tested to give correct values for j ≤ 128."
function Φ(n, q, j)
    if 0 ≤ n ≤ j
        return ϕ(n, q, j)
    elseif j < n ≤ 2 * j
        return (-1.0)^q * ϕ(2 * j - n, q, j)
    end
end

"`Φ2(n₁, n₂, q₁, q₂, j)` evaluates the two-dimensional Kravchuk functions over the set of lineal independent eigenmodes"
Φ2(n₁, n₂, q₁, q₂, j) = Φ(n₁, q₁, j) * Φ(n₂, q₂, j)
