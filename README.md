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
This file contains the function called ding which implements the projected conjugate gradient method proposed in the reference [2] to find zeros of a vector field. To invoke this function you need the following informations:

- x0 (vector{Float64}) current estimation to be a zero of the vector field;
- F (function) objective function.
- proj (function) maps a point x into its projection into a feasible set.

Moreover, you can specify the maximum number of iterates allowed, the tolerance given, and the value of the parameters xi, sigma, rho $\in (0,1)$, eta $\in [0,1)$ e theta $\in [0,1]$.

This function returns the following informations:

- x (vector{Float64}) current estimation to be a zero of the vector field;
- k (Int) number of iterations;
- t (Float64) elapsed time to solve the problem;
- nFx (Float64) contain the value of the euclidean norm of the current estimation;
- Fevals (Int) number of function evaluations;
- ierror (Int) the value stored in this variable tells the following messages: 0 - ||F(x0)||<tol, 1 - ||F(zk)||<tol, 2 - the maximum number of iterations has been exceeded.

# References 
[1] Surjanovic, S. & Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. Retrieved July 23, 2024, from http://www.sfu.ca/~ssurjano.  
[2] Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438.

