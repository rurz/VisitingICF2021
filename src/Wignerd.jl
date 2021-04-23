"This scripts computes the Wigner little-d functions for the rotations"

export dd, du

const β = π/2

"Lower triangle of modes. The counters are bounded as n ≤ 2j, and |m| ≤ n in steps of 2."

"coeff₁ to coeff₄ are the coefficients for the three-term recursions of the Wigner little-d"

coeff₁(n, k) = (2 * (-n/2 - (k - n/2 - 1) * cos(β)) * csc(β)) / √(k * (n - k + 1))

coeff₂(n, k) = √((k - 1) * (n + 2 - k) / (k * (n + 1 - k)))

coeff₃(n, k, l) = (2 * (1 - k + n/2 + (l - n/2 - 1) * cos(β)) * csc(β)) / √(l * (n - l + 1))

coeff₄(n, k) = √((k - 1) * (n + 2 - k) / (k * (n + 1 - k)))

"`dd(n, j)` is a function who depends on the energy lever _n_ and gives an _n × n_ little-d coefficient matrix. n ∈ [0, 2j] for j ≥ 0. Warning: j is tested to be accurate below 28, that is, matrices of 56 points are safe to be considered as accurate."
function dd(n, j)
    vdd = zeros((2 * j + 1, 2 * j + 1))
    vdd[1, 1] = cos(β/2)^n
    vdd[2, 1] = √n * (cos(β) - 1) * csc(β) * vdd[1, 1]
    for k in 2:n
        vdd[k + 1, 1] = coeff₁(n, k) * vdd[k, 1] - coeff₂(n ,k) * vdd[k - 1, 1]
    end
    for k in 1:(n + 1)
        for l in 1:n
            if l == 1
                vdd[k, l + 1] = coeff₃(n, k, l) * vdd[k, l]
            else
                vdd[k, l + 1] = coeff₃(n, k, l) * vdd[k, l] - coeff₄(n, l) * vdd[k, l - 1]
            end
        end
    end
    return vdd
end

# For j = 28 (N = 56), there's no overflow. Need to check

coeff₅(n, k, j) = (2 * ((1/2) * (n - 4 * j) - (k - 1 + (1/2) * (n - 4 * j)) * cos(β)) * csc(β)) / √(k * (4 * j + 1 - k - n))

coeff₆(n, k, j) = √(((k - 1) * (4 * j + 2 - k - n)) / (k * (4 * j + 1 - k - n)))

coeff₇(n, k, l, j) = (2 * (1 - k + (1/2) * (4 * j - n) + (l - 1 + (1/2) * (n - 4 * j)) * cos(β)) * csc(β)) / √(l * (4 * j + 1 - l - n))

coeff₈(n, k, j) = √(((k - 1) * (4 * j + 2 - k - n)) / (k * (4 * j + 1 - k - n)))

"`du(n, j)` is a function who depends on the energy lever _n_ and gives an _n × n_ little-d coefficient matrix. n ∈ [2j, 4j] for j ≥ 0. This modes are counted in reverse, n = 4j is the lowest, and n = 2j is the highest. Warning: j is tested to be accurate below 28, that is, matrices of 56 points are safe to be considered as accurate."
function du(n, j)
    vdu = zeros((2 * j + 1, 2 * j + 1))
    vdu[1, 1] = cos(β/2)^(4 * j - n)
    vdu[2, 1] = √(4 * j - n) * (cos(β) - 1) * csc(β) * vdu[1, 1]
    for k in 2:(4 * j - n)
        vdu[k + 1, 1] = coeff₅(n, k, j) * vdu[k, 1] - coeff₆(n, k, j) * vdu[k - 1, 1]
    end
    for k in 1:(4 * j + 1 - n)
        for l in 1:(4 * j - n)
            if l == 1
                vdu[k, l + 1] = coeff₇(n, k, l, j) * vdu[k, l]
            else
                vdu[k, l + 1] = coeff₇(n, k, l, j) * vdu[k, l] - coeff₈(n, l, j) * vdu[k, l - 1]
            end
        end
    end
    return vdu
end
