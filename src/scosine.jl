#   Another function with nontrivial groups and
#   repetitious elements.
#   NB: scaled version of COSINE

#   Source:
#   N. Gould, private communication.

#   classification OUR2-AN-V-0

#   Problem 50 in
#   L. Luksan, C. Matonoha and J. Vlcek
#   Modified CUTE problems for sparse unconstrained optimization,
#   Technical Report 1081,
#   Institute of Computer Science,
#   Academy of Science of the Czech Republic
#   http://www.cs.cas.cz/matonoha/download/V1081.pdf
#
# J.-P. Dussault, Clermont-Ferrand 05/2016.

# Note: discrepancy with CUTEst appears to be a bug in CUTEst, this matches the original paper
# (See issue #36)

export scosine, scosine_ADNLPModel

"Another function with nontrivial groups and repetitious elements in size 'n' "
function scosine(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("scosine: number of variables must be ≥ 2")
  n = max(2, n)

  p = zeros(n)
  for i=1:n
    p[i] = exp(6.0 * (i-1) / (n-1))
  end

  sum(cos(p[i]^2 * x[i]^2 - p[i+1] * x[i+1] / 2.0) for i=1:n-1)
end
function start_scosine(n :: Int)
  p = zeros(n)
  for i=1:n
    p[i] = exp(6.0 * (i-1) / (n-1))
  end
  x0 = map(pᵢ -> 1.0/pᵢ, p)
  return x0
end
scosine_ADNLPModel(n :: Int=100) = RADNLPModel(scosine, start_scosine(n), name="scosine "*string(n) * " variables")
