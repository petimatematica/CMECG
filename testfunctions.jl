# Reference: Surjanovic, S. & Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. Retrieved July 23, 2024, from http://www.sfu.ca/~ssurjano.

# Bowl-Shaped

# Sum Squares Function
function gradsumsquares(x)
    n = size(x)[1]
    G = Vector{Float64}(undef, n)
    for i in 1:n
        G[i]=2.0*i*x[i]
    end
    return G
end

function gradtrid(x) 
    n = size(x)[1]
    G = Vector{Float64}(undef, n)
    G[1]=2.0*(x[1]-1.0)-x[2]
    for i in 2:n-1
        G[i]=2.0*(x[i]-1.0)-x[i-1]-x[i+1]
    end
    G[n]=2.0*(x[n]-1.0)-x[n-1]
    return G
end

#=
Assumption 1 - 
Assumption 2 -
Assumption 3 - OK!=#
