#   Variable dimension diagonal quartic problem.
#
#   Source: problem 157 (p. 87) in
#   A.R. Buckley,
#   Test functions for unconstrained minimization,
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.
#
# See also
#
#   problem 23 in
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

export dqrtic, dqrtic_ADNLPModel

"Diagonal quartic model in size `n`"
function dqrtic(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)

  sum((x[i] - i)^4 for i=1:n)
end
start_dqrtic(n :: Int) =  (x -> 2 * x).(start_ones(n))
dqrtic_ADNLPModel(n :: Int=100) = RADNLPModel(dqrtic, start_dqrtic(n), name="dqrtic "*string(n) * " variables")
