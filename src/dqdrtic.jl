#   Diagonal quadratic problem
#
#   Source: problem 22 in
#   Ph. L. Toint,
#   "Test problems for partially separable optimization and results
#   for the routine PSPMIN",
#   Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
#
#   classification QUR2-AN-V-0

export dqdrtic, dqdrtic_ADNLPModel

"Diagonal quadratic problem"
function dqdrtic(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)

  sum(x[i]^2 + 100 * (x[i+1]^2 + x[i+2]^2) for i=1:n-2)
end
start_dqdrtic(n :: Int) =  (x -> 3 * x).(start_ones(n))
dqdrtic_ADNLPModel(n :: Int=100) = RADNLPModel(dqdrtic, start_dqdrtic(n), name="dqdrtic "*string(n) * " variables")
