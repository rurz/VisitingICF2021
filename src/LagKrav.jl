using VisitingICF2021

export  ϕd, ϕu, LKd, LKu

"`Φd(q₁, q₂, n, m, j)` shows the rotations of the Kravchuk function modes at (n, |m| ≤ n) for n ∈ [0, 2j]. It needs to be evaluted at (q₁, q₂) ∈ [0, 2j + 1]."
function Φd(q₁, q₂, n, m, j)
    v₀ = zeros(ComplexF64, 2 * j + 1)
    for n₁ in 0:n
        v₀[n₁ + 1] = Φ2(n₁, n - n₁, q₁, q₂, j) * dd(n, j)[n₁ + 1, Integer((m + n)/2) + 1]
    end
    return sum(v₀)
end

"`Φu(q₁, q₂, n, m, j)` shows the rotations of the Kravchuk function modes at (n, |m| ≤ 4j - n) for n ∈ [4j, 2j]. It needs to be evaluted at (q₁, q₂) ∈ [0, 2j + 1]."
function Φu(q₁, q₂, n, m, j)
    v₁ = zeros(ComplexF64, 2 * j + 1)
    for n₁ in 0:(4* j - n)
        v₁[n₁ + 1] = Φ2(2 * j - n₁, n + n₁ - 2 * j, q₁, q₂, j) * du(n, j)[n₁ + 1, Integer((m + 4 * j - n)/2) + 1]
    end
    return sum(v₁)
end

"Laguerre-Kravchuk modes functions are the core of the rotation of images, they encodes the original position information in a mode-angular momentum space, ready to apply an angle phase who rotates the image"

"`LKd(q₁, q₂, n, m, j)` shows the Laguerre-Kravchuk modes at (n, |m| ≤ n) for n ∈ [0, 2j]. It needs to be evaluted at (q₁, q₂) ∈ [0, 2j + 1]."
function LKd(q₁, q₂, n, m, j)
    vlkd = zeros(ComplexF64, 2 * j + 1)
    for n₁ in 0:n
        vlkd[n₁ + 1] = (-1.0im)^(n - n₁) * Φ2(n₁, n - n₁, q₁, q₂, j) * dd(n, j)[n₁ + 1, Integer((m + n)/2) + 1]
    end
    return (-1.0)^((1/2) * (abs(m) - m)) * sum(vlkd)
end

"`LKu(q₁, q₂, n, m, j)` shows the Laguerre-Kravchuk modes at (n, |m| ≤ 4j - n) for n ∈ [4j, 2j]. It needs to be evaluted at (q₁, q₂) ∈ [0, 2j + 1]."
function LKu(q₁, q₂, n, m, j)
    vlku= zeros(ComplexF64, 2 * j + 1)
    for n₁ in 0:(4* j - n)
        vlku[n₁ + 1] = (-1.0im)^(n - n₁) * Φ2(2 * j - n₁, n + n₁ - 2 * j, q₁, q₂, j) * du(n, j)[n₁ + 1, Integer((m + 4 * j - n)/2) + 1]
    end
    return (-1.0)^((1/2) * (abs(m) - m)) * sum(vlku)
end
