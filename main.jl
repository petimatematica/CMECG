using Distributions

include("search.jl"); include("testfunctions.jl"); include("signal.jl"); include("projections.jl")

#iter, t, x, Fx, Fnorm, erro = algorithm(ones(100), ARWHEAD, projRplus)

x = sparseSignalGenerator(n = 32, nze = 1) # x means x*
y = rand(Normal(), 8, 32)                  # y is a linear operator
A = y'*y
m = y*x
x0 = y'*m

iter, t, xx, Fx, mv, erro = signAlgorithm(x0, A, m, signalMonotoneNonlinearSystem,projRplus)