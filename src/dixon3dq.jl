#   Dixon's tridiagonal quadratic.
#
#   Source: problem 156 (p. 51) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.
#
#   classification QUR2-AN-V-0

export dixon3dq, dixon3dq_ADNLPModel

"Dixon's tridiagonal quadratic."
function dixon3dq(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  (x[1] - 1.0)^2 + (x[n] - 1.0)^2 + sum((x[i] - x[i+1])^2 for i=2:n-1)
end
start_dixon3dq(n :: Int) =  (xᵢ-> -1* xᵢ).(start_ones(n))
dixon3dq_ADNLPModel(n :: Int=10) = RADNLPModel(dixon3dq, start_dixon3dq(n), name="dixon3dq "*string(n) * " variables")
