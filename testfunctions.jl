# Reference: Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex
# constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438

using LinearAlgebra

function problemI(x)           # same name in the article

    v = Float64[]               # creates an empty Float64 vector

    for i in 1:length(x)        
        push!(v, exp(x[i])-1)   # calculates the value of each coordinate
    end

    return v                    # v means F(v) where F is a vector field

end

# problem_II
function problemII(x)                      # same name in the article

    v = Float64[]                           # creates an empty Float64 vector

    for i in 1:length(x)
        push!(v, x[i]-sin(abs(x[i]-1)))     # calculates the value of each coordinate
    end

    return v                                # v means F(v) where F is a vector field

end

# Remark: the DIAGONAL6 it is the problemI but unconstrained

function ARWHEAD(x)                                  # same name in the article
    
    v = Float64[]                                   # creates an empty Float64 vector
    n = length(x)

    for i in 1:n-1
        push!(v, -4+4*x[i]*(x[i]^2+x[n]^2)) # calculates the value of each coordinate, except the last
    end

    push!(v, 4*x[n]*sum(x[j]^2+x[n]^2 for j in 1:n-1)) # calculates the last coordinate

    return v                                        # v means F(v) where F is a vector field

end

function PENALTY1(x; c=1.e-5)                      # same name in the article

    v = Float64[]                           # creates an empty Float64 vector
    n = length(x)

    for i in 1:n
        push!(v, 2*c*(x[i]-1)+4*(sum(x[j]^2 for j in i:n)-0.25)*x[i])
    end

    return v                                # v means F(v) where F is a vector field

end