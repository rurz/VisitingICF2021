using VisitingICF2021
using LinearAlgebra

function fnm(n, m, j, imgdat)
    vnm = zeros(ComplexF64, (2 * j + 1, 2 * j + 1))
    for q₁ in -j:j
        for q₂ in -j:j
            vnm[q₁ + j + 1, q₂ + j + 1] = LK(q₁, q₂, n, m, j) * imgdat[q₁ + j + 1, q₂ + j + 1]
        end
    end
    return sum(vnm)
end

fnmθ(θ, n, m, j, imgdat) = exp(-1im * m * θ) * fnm(n, m, j, imgdat)

function Fqθ(θ, q₁, q₂, j, imgdat)
    rd = 0.0
    ru = 0.0
    rd = sum([conj(LK(q₁, q₂, n, m, j)) * fnmθ(θ, n, m, j, imgdat) for n in 0:2*j for m in -n:2:n])
    ru = sum([conj(LK(q₁, q₂, n, m, j)) * fnmθ(θ, n, m, j, imgdat) for n in 2*j:4*j for m in -(4 * j - n):2:(4 * j - n)])
    return rd + ru
end
