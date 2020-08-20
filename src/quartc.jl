#   A simple quartic function.
#
#   Source:  problem 157 (p. 87) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.
#
#   classification OUR2-AN-V-0

export quartc, quartc_ADNLPModel

"A simple quartic function."
function quartc(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)

  sum((x[i] - i)^4 for i=1:n)
end
start_quartc(n :: Int) =  (x -> 2 * x).(start_ones(n))
quartc_ADNLPModel(n :: Int=100) = RADNLPModel(quartc, start_quartc(n), name="quartc "*string(n) * " variables")
