include("procedures.jl"); include("testfunctions.jl"); include("projections.jl")
using CUTEst, NLPModels
#A, B, C, p, q, r = myresults(gradsphere, id, factor =50)

#nlp = CUTEstModel("PENALTY1","-param","N=100")


x0 = -1.0*ones(90)

function F(x)
    #return grad(nlp,x)
    
    n = length(x)
    G = Vector{Float64}(undef,n)

    for i in 1 : n
        G[i] = 2. * (n - i + 1.) * x[i]
    end

    return G

end

s,iters,et,norm_F_x, error = ding(x0, F, id)


finalize(nlp)