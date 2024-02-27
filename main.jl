using Distributions, Random

include("search.jl"); include("testfunctions.jl"); include("projections.jl")

n = 4096

# Dimensão do sinal eh 4n

global A = rand(Normal(),Int(n/4),n)

rng = MersenneTwister(1234)

sinal_original = rand(rng,n)

global b = A * sinal_original

y = A' * b

e2n = ones(2n)

c = 0.1 * e2n + [-y;y]

ATA = A' * A
H = [ATA -ATA ; -ATA ATA]

function F(z)

    d = H * z + c
    return min.(z,d)
end

x0 = ones(2n) 


iter, t, x, Fx, Fnorm, erro = algorithm(x0,F, projRplus)

#x = sparseSignalGenerator(n = 32, nze = 1) # x means x*
#y = rand(Normal(), 8, 32)                  # y is a linear operator
#A = y'*y
#m = y*x
#x0 = y'*m
#
#iter, t, xx, Fx, mv, erro = signAlgorithm(x0, A, m, signalMonotoneNonlinearSystem,projRplus)



