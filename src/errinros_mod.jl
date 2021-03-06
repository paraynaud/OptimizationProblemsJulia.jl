# Errin Rosenbrock - modified function.

#   problem 27 in
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
# J,-P, Dussault, Rennes 09/2015.

export errinros_mod, errinros_mod_ADNLPModel


function errinros_mod(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)
  n < 2 && @warn("errinros_mod: number of variables must be ≥ 2")
  n = max(2, n)

  sum((x[i-1] - 16.0 * x[i]^2 * (1.5 + sin(i))^2)^2 for i=2:n) + sum((1.0 - x[i])^2 for i=2:n)
end
start_errinros_mod(n :: Int) =  (x -> -1 * x).(start_ones(n))
errinros_mod_ADNLPModel(n :: Int=100) = RADNLPModel(errinros_mod, start_errinros_mod(n), name="errinros_mod "*string(n) * " variables")
