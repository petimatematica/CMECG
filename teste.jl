using CUTEst, NLPModels

include("search.jl"); include("testfunctions.jl"); include("projections.jl")

problems = CUTEst.select()

problem = problems[16]

nlp = CUTEstModel(problem)

x0 = nlp.meta.x0 # armazena o chute inicial sugerido na variável x0

g(x)=grad(nlp, x) # definindo o gradiente

iter, t, x, Fx, Fnorm = algorithm(x0, g, id)

finalize(nlp)

## Solved problems CUTEst ##

# 8 - LUKSAN13LS; 11 - QPCBOEI1