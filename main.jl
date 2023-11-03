include("search.jl"); include("testfunctions.jl")

iter, t, x, Fx, Fnorm, erro = algorithm(ones(100), gradsphere, id)