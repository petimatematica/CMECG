include("procedures.jl"); include("testfunctions.jl"); include("projections.jl")

x0 = ones(100) # guess
F = gradsumsquares # vector field
proj = projRplus # projection

# x,k,t,nFx,Fevals,erro=ding(x0,F,proj, maxiter=1.e4) 
# use the line above if you want to solve the equation F(x)=0 for more details see Readme.md
