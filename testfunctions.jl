# Reference: Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex
# constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438

using LinearAlgebra

# problem_I
function problem_I(x)
    v = Float64[]
    for i in 1:length(x)
        push!(v, exp(x[i])-1)
    end

    return v

end

function projproblem_I(x)

    v = Float64[]

    for i in 1:length(x)
        if x[i] < 0
            push!(v,0.0)
        else
            push!(v, x[i])
        end
    end

    return v

end

function projhyperplane(y; a=ones(length(y)))
    return y + (length(y)-dot(a,y))/norm(a)^2*a
end

function ncheck(x)
    flag = 0
    for i in 1:length(x)
        if x[i] < 0
            return 1
        end
    end
    return flag
end

function vlist(n)
    A = Array{Vector}(undef, n+1)
    A[1] = zeros(n)
    for i in 1:n
        v = zeros(n)
        v[i] = n
        A[i+1] = v
    end
    return A
end

function leastdist(y)
    x = Array{Vector}(undef, 3)
    x[1] = y
    n = length(y)
    A = vlist(n)

    if norm(y-A[1]) < norm(y-A[2])
        x[2] = A[1]
        x[3] = A[2]
    else
        x[2] = A[2]
        x[3] = A[1]
    end

    for i in 3:length(A)
        if norm(y-A[i]) < norm(y-x[2])
            x[3] = x[2]
            x[2] = A[i]
        elseif norm(y-A[i]) < norm(y-x[3]) && norm(y-A[i]) > norm(y-x[2])
            x[3] = A[i]
        end
    end

    return x
end


# problem_II
function problem_II(x)
    v = Float64[]
    for i in 1:length(x)
        push!(v, x[i]-sin(abs(x[i]-1)))
    end

    return v

end

function projproblem_II(x)
    n = length(x)
    y = projproblem_I(x)
    if sum(y) > n || ncheck(y) == 1
        y = projhyperplane(y,ones(n))
        if ncheck(y) == 1
            A = leastdist(y)
            u = y-A[2]
            v = A[3]-A[2]
            if dot(u,v)/(norm(u)*norm(v)) < 0
                y = A[2]
            else
                y = dot(v,u)/dot(v,v)*v
            end
        end
    end 
    return y
end

id(x) = x

function gradsphere(x)
    v = Float64[]
    for i in 1:length(x)
        push!(v, 2x[i])
    end

    return v
    
end

function projsphere(y; x = zeros(length(y)), δ = 1)
    d = norm(x-y)
    if d > δ
        return x + δ*(y-x)/d
    else
        return y
    end
end

# DIAGONAL6
function diagonal6(x)
    v = Float64[]
    for i in 1:length(x)
        push!(v, exp(x[i])-1)
    end

    return v

end

# ARWHEAD
function arwhead(x)
    v = Float64[]
    l = length(x)
    for i in 1:l
        if i == l
            push!(v, 4*x[i]*sum(x[j]^2+x[i]^2 for j in 1:(l-1)))
        else
            push!(v, -4+4*x[i]*(x[i]^2+x[l]^2))
        end
    end
    return v
end

# PENALTYI
function penalty_I(x; c = 1.0e-5)
    v = Float64[]
    n = length(x)
    for i in 1:n
        push!(v, 2*c*(x[i]-1)+4*(sum(x[j]^2 for j in 1:n)-0.25)*x[i])
    end
    return v
end

# DIXON3DQ
function dixon3dq(x)
    v = Float64[]
    n = length(x)
    m = Int(floor((n-1)/2))
    push!(v, 2*(x[1]-1))

    for i in 2:(n-1)
        if i%2 == 0
            push!(v, 2*sum(x[k]-x[k+1] for k in 2:2*m))
        else
            push!(v, -2*sum(x[k]-x[k+1] for k in 2:2*m+1))
        end
    end

    push!(v, 2*(x[n]-1))

    return v

end

# GENHUMPS
function genhumps(x)
    v = Float64[]
    n = length(x)
    m = Int(floor((n-1)/2))

    for i in 1:n
        if i%2 == 0
            push!(v, 2*sum(sin(4*x[k]*sin(2*x[k-1])^2) for k in 2:2*m) + 0.1*sum(x[k] for k in 2:2*m))
        else
            push!(v, 2*sum(sin(4*x[k]*sin(2*x[k-1])^2) for k in 2:2*m+1) + 0.1*sum(x[k] for k in 2:2*m+1))
        end
    end

    return v

end

# ENGVALI
function engvali(x)

    v = Float64[]
    n = length(x)
    push!(v, 4*x[1]*(x[1]^2+x[2]^2)-4)

    for i in 2:n-2
        push!(v, 4*x[i]*(x[i-1]^2+x[i]^2)+4*x[i]*(x[i]^2+x[i+1]^2)-4)
    end

    push!(v, 4*x[n-1]*(x[n-2]^2+x[n-1]^2))
    push!(v, 4*x[n]*(x[n-1]^2+x[n]^2))

    return v

end

# DIXMAANH
function dixmaanh(x; α = 1, β = 0.26, γ = 0.26, δ = 0.26, k1 = 1, k4 = 1, k2 = 0, k3 = 0)

    v = Float64[]
    n = length(x)
    m = Int(floor(n/3))
    push!(v, 2*α*x[1]*(1/n)^(k1) + 2*β*x[1]*(x[2]+x[2]^2)^2*(1/n)^(k2) + 2*γ*x[1]*x[1+m]^4*(1/n)^(k3) + δ*x[1+2m]*(1/n)^(k4x))
    
    for i in 2:m
        push!(v, 2*α*x[i]*(i/n)^(k1) + 2*β*(x[i-1]^2)*(x[i]+x[i]^2)^2*(1+2*x[i])*((i-1)/n)^(k2) + 2*β*x[i]*(x[i+1]+x[i+1]^2)^2*(i/n)^(k2) + 2*γ*x[i]*x[i+m]^4*(1/n)^(k3) + δ*x[i+2m]*(1/n)^(k4))
    end

    for i in m+1:2*m
        push!(v, 2*α*x[i]*(i/n)^(k1) + 2*β*(x[i-1]^2)*(x[i]+x[i]^2)^2*(1+2*x[i])*((i-1)/n)^(k2) + 2*β*x[i]*(x[i+1]+x[i+1]^2)^2*(i/n)^(k2) + 2*γ*x[i]*x[i+m]^4*(1/n)^(k3) + 4*γ*x[i-m]^2*x[i]^3*((i-m)/n)^(k3))
    end

    for i in 2*m:n-1
        push!(v, 2*α*x[i]*(i/n)^(k1) + 2*β*(x[i-1]^2)*(x[i]+x[i]^2)^2*(1+2*x[i])*((i-1)/n)^(k2) + 2*β*x[i]*(x[i+1]+x[i+1]^2)^2*(i/n)^(k2) + 4*γ*x[i-m]^2*x[i]^3*((i-m)/n)^(k3) + δ*x[i-2m]*((i-2m)/n)^(k4))
    end

    push!(v, 2*α*x[n])
end