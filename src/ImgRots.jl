"This script contains the implementation algorithms for the unitary transformation on discrete and finite fields. It acts as a single transformation, but it can be adecuate to be extended on a FOR loop to be applied on polychromatic fields

It provides three functions: `fnm(n, m, j, imgdata)` maps the field from position to modes. `fnmθ(θ, n, m, j, imgdat)` adds an angle dependent phase. `Fqθ(θ, q₁, q₂, j, imgdat)` returns the transformed field in the modes to the position representation. This last function is the upper layer of the transformation, it can (could, should) be used to the polychromatic duties."

using VisitingICF2021 # Call the parent module

export fnm, fnmθ, Fqθ

"`fnm(n, m, j, imgdat)` encodes a two-dimensional position field in a modes space labeled by mode _n_ and counter _m_, who ranges are n ∈ [0, 2j] and |m|₂ ≤ n. 'imgdat' is a position field, that is a N×N matrix, for N = 2j."
function fnm(n, m, j, imgdat)
    vnm = zeros(ComplexF64, (2 * j + 1, 2 * j + 1))
    for q₁ in -j:j
        for q₂ in -j:j
            vnm[q₁ + j + 1, q₂ + j + 1] = LK(q₁, q₂, n, m, j) * imgdat[q₁ + j + 1, q₂ + j + 1]
        end
    end
    return sum(vnm)
end

"`fnmθ(θ, n, m, j, imgdat)` is a function who applies an angle-dependent phase to the `fnm` function. It is the preamble to perform the transformation, given the value of θ."
fnmθ(θ, n, m, j, imgdat) = exp(-1im * m * θ) * fnm(n, m, j, imgdat)

"`Fqθ(θ, q₁, q₂, j, imgdat)` is the top function who perform the whole transformation on the position field. Given and N×N field, an angle θ, the values of the transformed pair (q₁, q₂) are given by this function. Note: It is extremely slow, so slow."
function Fqθ(θ, q₁, q₂, j, imgdat)
    rd = 0.0
    ru = 0.0
    rd = sum([conj(LK(q₁, q₂, n, m, j)) * fnmθ(θ, n, m, j, imgdat) for n in 0:2*j for m in -n:2:n])
    ru = sum([conj(LK(q₁, q₂, n, m, j)) * fnmθ(θ, n, m, j, imgdat) for n in (2*j + 1):4*j for m in -(4 * j - n):2:(4 * j - n)])
    return rd + ru
end
