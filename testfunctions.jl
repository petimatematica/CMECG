#= Reference: Surjanovic, S. & Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. 
              Retrieved July 23, 2024, from http://www.sfu.ca/~ssurjano.=#

# Bowl-Shaped

# Rotated Hyper-Ellipsoid Function
function gradrotated(x)
    n = size(x)[1]
    G = Vector{Float64}(undef, n)
    for i in 1:n
        G[i]=2.0*(n-1+1.0)*x[i]
    end
    return G
end

# Sphere Function
gradsphere(x)=2x

# Sum Squares Function
function gradsumsquares(x)
    n = size(x)[1]
    G = Vector{Float64}(undef, n)
    for i in 1:n
        G[i]=2.0*i*x[i]
    end
    return G
end
