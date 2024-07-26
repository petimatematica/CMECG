include("testfunctions.jl"); include("projections.jl")

x0 = ones(100) # guess
F = gradsumsquares # see testfunctions.jl for more details
proj = projRplus # see projections.jl for more details

x,k,t,nFx,Fevals,error=ding(x0,F,proj, maxiter=1.e4) # see search.jl for more details
