include("search.jl"); include("testfunctions.jl"); include("projections.jl")

x0 = ones(100)

iter, t, x, Fx, normFx, stpcond = algorithm(x0, ackley, projRplus) # try to find the x in n dimensial euclidean space, such that F(x) = 0.

# iter: number of iterates
# t: elapsed time
# x: an element of the n dimensional euclidean space
# Fx: the value of the objective function at the point x
# normFx: norm of Fx (Here is where we can check if we find a solution)
# stpcond: indicates which stop condition holds. If stpcond = 0, then ||f(x)|| < \epsilon and if stpcond = 1, we have ||f(z)|| < \epsilon
