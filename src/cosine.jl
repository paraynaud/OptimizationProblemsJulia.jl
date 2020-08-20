# The cosine function.
#
#   Source: problem 6 in
#   L. Luksan, C. Matonoha and J. Vlcek
#   Modified CUTE problems for sparse unconstrained optimization,
#   Technical Report 1081,
#   Institute of Computer Science,
#   Academy of Science of the Czech Republic
#
#   http://www.cs.cas.cz/matonoha/download/V1081.pdf
#
#   classification OUR2-AN-V-0
#
# D. Orban, Montreal, 08/2015.

export cosine, cosine_ADNLPModel

"The cosine function in size `n`"
function cosine(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("cosine: number of variables must be ≥ 2")
  n = max(2, n)
  sum(cos(x[i]^2 - 0.5 * x[i+1]) for i = 1:n-1)
end
start_cosine(n :: Int) = start_ones(n)
cosine_ADNLPModel(n :: Int=100) = RADNLPModel(cosine, start_cosine(n), name="cosine "*string(n) * " variables")
