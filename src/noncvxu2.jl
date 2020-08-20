#   A nonconvex unconstrained function with a unique minimum value

#   classification OUR2-AN-V-0

#   Problem 43 in
#   L. Luksan, C. Matonoha and J. Vlcek
#   Modified CUTE problems for sparse unconstrained optimization,
#   Technical Report 1081,
#   Institute of Computer Science,
#   Academy of Science of the Czech Republic
#   http://www.cs.cas.cz/matonoha/download/V1081.pdf
#
# J.-P. Dussault, Clermont-Ferrand 05/2016.

export noncvxu2, noncvxu2_ADNLPModel

function noncvxu2(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("noncvxu2: number of variables must be ≥ 2")
  n = max(2, n)

  sum((x[i] + x[mod(3 * i - 2, n) + 1] + x[mod(7 * i - 3, n) + 1])^2 +
  4.0 * cos(x[i] + x[mod(3 * i - 2, n) + 1] + x[mod(7 * i - 3, n) + 1]) for i=1:n)
end
start_noncvxu2(n :: Int) =  (x -> (Float64)(x)).([1:n;])
noncvxu2_ADNLPModel(n :: Int=100) = RADNLPModel(noncvxu2, start_noncvxu2(n), name="noncvxu2 "*string(n) * " variables")
