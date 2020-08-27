# A banded function with semi-bandwidth b and
# negative curvature near the starting point.
#
# Note that the initial point in the reference below is erroneous.
# In this model, we use the starting point specified in the
# original SIF model, part of the CUTE collection.
#
# See also
#
#   problems 8, 9, 10 in
#   L. Luksan, C. Matonoha and J. Vlcek
#   Modified CUTE problems for sparse unconstrained optimization,
#   Technical Report 1081,
#   Institute of Computer Science,
#   Academy of Science of the Czech Republic
#
#   http://www.cs.cas.cz/matonoha/download/V1081.pdf
#
#   classification SUR2-AN-V-0
#
# D. Orban, Montreal, 08/2015.

export curly, curly10, curly20, curly30, curly_ADNLPModel, curly10_ADNLPModel, curly20_ADNLPModel, curly30_ADNLPModel

curly10(x :: AbstractVector{Y}) where Y <: Number = curly(x, b=10)
curly20(x :: AbstractVector{Y}) where Y <: Number = curly(x, b=20)
curly30(x :: AbstractVector{Y}) where Y <: Number = curly(x, b=30)

"Curly function in size `n` with semi-bandwidth `b`"
function curly(x :: AbstractVector{Y}; b :: Int=10) where Y <: Number
  n = length(x)
  n < 2 && @warn("curly: number of variables must be ≥ 2")
  n = max(2, n)

  f = Vector{Y}(undef,n)
  map!(i -> sum(x[j] for j=i:min(i+b,n)),f, [1:n;])

  sum(f[i] * (f[i] * (f[i]^2 - 20) - 0.1) for i = 1:n)
end
start_curly(n :: Int) = [1.0e-4 * i /(n+1) for i = 1:n]
curly_ADNLPModel(n :: Int=100) = RADNLPModel(curly, start_curly(n), name="curly "*string(n) * " variables")
curly10_ADNLPModel(n :: Int=100) = RADNLPModel(curly10, start_curly(n), name="curly10 "*string(n) * " variables")
curly20_ADNLPModel(n :: Int=100) = RADNLPModel(curly20, start_curly(n), name="curly20 "*string(n) * " variables")
curly30_ADNLPModel(n :: Int=100) = RADNLPModel(curly30, start_curly(n), name="curly30 "*string(n) * " variables")
