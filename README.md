# main.jl
This file contains the main file which includes all the other files presented in this readme. Its code is written below.

```julia
include("testfunctions.jl"); include("projections.jl")

x0 = ones(100) # guess
F = gradsumsquares # see testfunctions.jl for more details
proj = projRplus # see projections.jl for more details

x,k,t,nFx,Fevals,erro=ding(x0,F,proj, maxiter=1.e4) # see search.jl for more details
```
# testfunction.jl
This file contain the test functions used in our work.

# testfunction.jl
This file contain the projections used in our work.

# search.jl
This file contains the function called ding which implements the PRP method (conjugate gradient) or the gradient method where the steplength is computated by Armijo's rule. To invoke this function you need the following informations:

# References 
[1] Surjanovic, S. & Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. Retrieved July 23, 2024, from http://www.sfu.ca/~ssurjano.  
[2] Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438

