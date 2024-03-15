# Reference:
# Surjanovic, S. & Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. 
# Retrieved March 15, 2024, from http://www.sfu.ca/~ssurjano.

using LinearAlgebra
# ackley function

function ackley(x; a=20, b=0.2, c=2pi)
    
    d = length(x)
    
    x1 = -a*exp(-b*sqrt(1/d*sum(x[i]^2 for i in 1:d)))
    x2 = -exp(1/d*sum(cos(c*x[i]) for i in 1:d))
    
    return x1 + x2 + a + exp(1)

end

function gradackley(x; a=20, b=0.2, c=2pi)

    v = Float64[]

    d = length(x)
    k = sqrt(1/d*sum(x[i]^2 for i in 1:d))
    
    for i in 1:d
        push!(v, -(a*b)/(d*k)*exp(-b*k)*x[i]+exp(1/d*sum(cos(c*x[j]) for j in 1:d)) + c/d*sin(c*x[i]))
    end

    return v

end

# trid function

trid(x) = sum((x[i]-1)^2 for i in 1:size(x)[1]) - sum(x[j]*x[j-1] for j in 2:size(x)[1])

function gradtrid(x)
    
    d = length(x)
    v = Float64[]

    push!(v,2*(x[1]-1)-x[2])

    for i in 2:d-1
        push!(v, 2*(x[i]-1)-x[i-1]-x[i+1])
    end

    push!(v, 2*(x[d]-1)-x[d-1])
    
    return v

end

# zakharov function
function zak(x)
    d = size(x)[1]
    x1 = sum(x[i]^2 for i in 1:d)
    x2 = sum(0.5*j*x[j] for j in 1:d)^2

    return x1 + x2 + x2^2
end

function gradzak(x)
    d = size(x)[1]
    k = sum(0.5*j*x[j] for j in 1:d)
    v = Float64[]
    for i in 1:d
        push!(v, 2*x[i] + i*k + 2*i*k^3)
    end
    return v
end

# Dixon-Price Function

dixon(x) = (x[1]-1)^2 + sum(i*(2*x[i]^2-x[i-1])^2 for i in 2:length(x))

function graddixon(x)
    d = size(x)[1]
    v = Float64[]
    push!(v, 2*(x[1]-1)-4*(2*x[2]^2-x[1]))
    for i in 2:d-1
        push!(v, 8*i*x[i]*(2*x[i]^2-x[i-1])-2*(i+1)*(2*x[i+1]^2-x[i]))
    end
    push!(v, 8*d*x[d]*(2*x[d]^2-x[d-1]))
    return v
end

# Rosembrock function

rosembrock(x) = sum(100*(x[i+1]-x[i]^2)^2+(x[i]-1)^2 for i in 1:length(x)-1)

function gradrosembrock(x)
    d = size(x)[1]
    v = Float64[]
    push!(v, 2*(x[1]-1)-400*x[1]*(x[2]-x[1]^2))
    for i in 2:d-1
        push!(v, 200*(x[i]-x[i-1]^2)-400*x[i]*(x[i+1]-x[i]^2)+2*(x[i]-1))
    end
    push!(v, 200*(x[d]-x[d-1]^2))
    return v
end
