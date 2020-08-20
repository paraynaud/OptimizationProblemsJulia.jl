using Test, OptimizationProblemsJulia


# Test that every problem can be instantiated.
all_probs = names(OptimizationProblemsJulia)
f(x :: Symbol) = begin string_x = string(x); length(string_x) >=10 ? string_x[end-9:end] == "ADNLPModel" : false end
probs = filter(x-> f(x), all_probs)
for prob in probs
  prob == :OptimizationProblemsJulia && continue
  println(prob)
  # @show typeof(prob)
  prob_fn = eval(prob)
  prob_fn()
end
