#   Source:
#   Ph. Toint, private communication, 1997.

#   SDIF input: N. Gould and Ph. Toint, November 1997.

#   problem 35 in
#   L. Luksan, C. Matonoha and J. Vlcek
#   Modified CUTE problems for sparse unconstrained optimization,
#   Technical Report 1081,
#   Institute of Computer Science,
#   Academy of Science of the Czech Republic
#
#   http://www.cs.cas.cz/matonoha/download/V1081.pdf
#

#   classification SUR2-AN-V-0
# J.-P. Dussault, Rennes 09/2015.

export genhumps, genhumps_ADNLPModel

function genhumps(x :: AbstractVector{Y}) where Y <: Number
  n = length(x)

  ζ = 20.0
  sum(( sin(ζ * x[i])^2 * sin(ζ * x[i+1])^2 + 0.05 * (x[i]^2 + x[i+1]^2)) for i=1:n-1)
end
start_genhumps(n :: Int) =  begin x0 = (x -> -506.2 * x).(start_ones(n)); x0[1] = - 506.0; return x0 end
genhumps_ADNLPModel(n :: Int=100) = RADNLPModel(genhumps, start_genhumps(n), name="genhumps "*string(n) * " variables")
