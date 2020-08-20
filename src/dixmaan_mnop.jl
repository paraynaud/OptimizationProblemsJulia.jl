
# The Dixon-Maany test problem (version M by default)
#
#   Source:
#   L. C. W. Dixon and Z. Maany,
#   A family of test problems with sparse Hessians for unconstrained
#   optimization,
#   TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.
#
# See also
#
#   problems 19, 20, 21, 22 in
#   L. Luksan, C. Matonoha and J. Vlcek
#   Modified CUTE problems for sparse unconstrained optimization,
#   Technical Report 1081,
#   Institute of Computer Science,
#   Academy of Science of the Czech Republic
#
#   http://www.cs.cas.cz/matonoha/download/V1081.pdf
#
# D. Orban, Montreal, 08/2015.

export dixmaanm, dixmaann, dixmaano, dixmaanp, dixmaanm_ADNLPModel, dixmaann_ADNLPModel, dixmaano_ADNLPModel, dixmaanp_ADNLPModel

dixmaann(x :: AbstractVector{Y}) where Y <: Number = dixmaanm(x, α=1.0, β=0.0625, γ=0.0625, δ=0.0625)
dixmaano(x :: AbstractVector{Y}) where Y <: Number = dixmaanm(x, α=1.0, β=0.125, γ=0.125, δ=0.125)
dixmaanp(x :: AbstractVector{Y}) where Y <: Number = dixmaanm(x, α=1.0, β=0.26, γ=0.26, δ=0.26)

"Dixon-Maany function in size `n` (version M by default)"
function dixmaanm(x :: AbstractVector{Y};
                  α :: Float64=1.0, β :: Float64=0.0, γ :: Float64=0.125, δ :: Float64=0.125) where Y <: Number
  n = length(x)
  (n % 3 == 0) || @warn("dixmaan: number of variables adjusted to be a multiple of 3")
  m = max(1, div(n, 3))
  n = 3 * m

  1 +
  sum((i / n)^2 * α * x[i]^2 for                     i=1:n) +
  sum(i / n * β * x[i]^2 * (x[i+1] + x[i+1]^2)^2 for i=1:n-1) +
  sum(i / n * γ * x[i]^2 * x[i+m]^4 for              i=1:2*m) +
  sum((i / n)^2 * δ * x[i] * x[i+2*m] for            i=1:m)

end
start_dixmaanm(n :: Int) = [2.0 for i = 1:n]
dixmaanm_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanm, start_dixmaanm(n), name="dixmaanm "*string(n) * " variables")
dixmaann_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaann, start_dixmaanm(n), name="dixmaann "*string(n) * " variables")
dixmaano_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaano, start_dixmaanm(n), name="dixmaano "*string(n) * " variables")
dixmaanp_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanp, start_dixmaanm(n), name="dixmaanp "*string(n) * " variables")
