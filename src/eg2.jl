# A simple non convex problem with several local minima.
#
#   Source: Section 1.2.4 of
#   A. R. Conn, N. I. M. Gould and Ph. L. Toint,
#   LANCELOT, A Fortran Package for Large-Scale Nonlinear Optimization
#   (Release A)
#   Springer Verlag, 1992.
#
# See also
#
#   problem 25 in
#   L. Luksan, C. Matonoha and J. Vlcek
#   Modified CUTE problems for sparse unconstrained optimization,
#   Technical Report 1081,
#   Institute of Computer Science,
#   Academy of Science of the Czech Republic
#
#   http://www.cs.cas.cz/matonoha/download/V1081.pdf
#
#   classification OUR2-AN-1000-0
#
# D. Orban, Montreal, 08/2015.

export eg2, eg2_ADNLPModel

"model in size `n`"
function eg2(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("eg2: number of variables must be ≥ 2")
  n = max(2, n)

    sum(
      sin(x[1] + x[i]^2 - 1)
      for i=1:n-1
    ) +
    0.5 * sin(x[n]^2)
end
start_eg2(n :: Int) =  (x -> 0 * x).(start_ones(n))
eg2_ADNLPModel(n :: Int=100) = RADNLPModel(eg2, start_eg2(n), name="eg2 "*string(n) * " variables")
