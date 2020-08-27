using Test, OptimizationProblemsJulia
using OptimizationProblems
using NLPModels, NLPModelsJuMP
using LinearAlgebra

# Test that every problem can be instantiated.
all_probs = names(OptimizationProblemsJulia)
f(x :: Symbol) = begin string_x = string(x); length(string_x) >=10 ? string_x[end-9:end] == "ADNLPModel" : false end
probs = filter(x-> f(x), all_probs)
f_jump(x :: Symbol) = begin string_x = string(x); return Symbol(string_x[1:end-11]) end

σ(a,b) = abs(a-b) < 1e-6
τ(a,b) = norm(a-b) < 1e-6


@testset "test on default length/initial point/obj/grad/hv" begin
  for prob in probs
    prob == :OptimizationProblemsJulia && continue
    prob_jump = f_jump(prob)

    prob_fn = eval(prob)
    prob_fn_jump = OptimizationProblems.eval(prob_jump)

    ad_nlp = prob_fn()
    jump_nlp = NLPModelsJuMP.MathOptNLPModel(prob_fn_jump())

    n_ad = ad_nlp.meta.nvar
    n_jump = jump_nlp.meta.nvar

    ad_x0 = ad_nlp.meta.x0
    jump_x0 = jump_nlp.meta.x0

    @test n_ad == n_jump
    @test ad_x0 == jump_x0

    x = ad_x0
    v = ones(n_ad)

    @test σ(NLPModels.obj(ad_nlp, x), NLPModels.obj(jump_nlp, x))
    @test τ(NLPModels.grad(ad_nlp, x), NLPModels.grad(jump_nlp, x))
    @test τ(NLPModels.hprod(ad_nlp, x, v), NLPModels.hprod(jump_nlp, x, v))

  end
end
