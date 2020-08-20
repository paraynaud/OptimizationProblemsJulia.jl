# The Dixon-Maany test problem (version I by default)
#
#   Source:
#   L. C. W. Dixon and Z. Maany,
#   A family of test problems with sparse Hessians for unconstrained
#   optimization,
#   TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.
#
# See also
#
#   problems 15, 16, 17, 18 in
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

export dixmaani, dixmaanj, dixmaank, dixmaanl, dixmaani_ADNLPModel, dixmaanj_ADNLPModel, dixmaank_ADNLPModel, dixmaanl_ADNLPModel


dixmaanj(x :: AbstractVector{Y}) where Y <: Number = dixmaani(x, α=1.0, β=0.0625, γ=0.0625, δ=0.0625)
dixmaank(x :: AbstractVector{Y}) where Y <: Number = dixmaani(x, α=1.0, β=0.125, γ=0.125, δ=0.125)
dixmaanl(x :: AbstractVector{Y}) where Y <: Number = dixmaani(x, α=1.0, β=0.26, γ=0.26, δ=0.26)

"Dixon-Maany function in size `n` (version I by default)"
function dixmaani(x :: AbstractVector{Y};
                  α :: Float64=1.0, β :: Float64=0.0, γ :: Float64=0.125, δ :: Float64=0.125) where Y <: Number
  n = length(x)
  (n % 3 == 0) || @warn("dixmaan: number of variables adjusted to be a multiple of 3")
  m = max(1, div(n, 3))
  n = 3 * m

  1 +
  sum((i / n)^2 * α * x[i]^2 for             i=1:n)   +
  sum(β * x[i]^2 * (x[i+1] + x[i+1]^2)^2 for i=1:n-1) +
  sum(γ * x[i]^2 * x[i+m]^4 for              i=1:2*m) +
  sum((i / n)^2 * δ * x[i] * x[i+2*m] for    i=1:m)
end
start_dixmaani(n :: Int) = [2.0 for i = 1:n]
dixmaani_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaani, start_dixmaani(n), name="dixmaani "*string(n) * " variables")
dixmaanj_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanj, start_dixmaani(n), name="dixmaanj "*string(n) * " variables")
dixmaank_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaank, start_dixmaani(n), name="dixmaank "*string(n) * " variables")
dixmaanl_ADNLPModel(n :: Int=99) = RADNLPModel(dixmaanl, start_dixmaani(n), name="dixmaanl "*string(n) * " variables")
