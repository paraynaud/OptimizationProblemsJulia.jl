module OptimizationProblemsJulia

# using OptimizationProblems
using ADNLPModels


start_ones(n :: Int) = ones(n)


path = dirname(@__FILE__)
files = filter(x->x[end-2:end] == ".jl", readdir(path))
for file in files
  @show file
  if (file == "OptimizationProblemsJulia.jl") || (file == "template.jl"); continue; end
  include(file)
end

#=
  à Faire :
  - schmvett
  - tquartic
=#

end # module
