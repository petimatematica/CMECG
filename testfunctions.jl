# Reference: https://www.sfu.ca/~ssurjano/optimization.html

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
