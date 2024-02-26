using LinearAlgebra, Distributions, SparseArrays

# Function: sparseSignalGenerator
# Description: Generate a sparse matrix using the given values of number of rows, columns
# and non-zero entries, respectively.

#----------Inputs---------------
# k - number of rows (default: 4096)
# n - number of columns (default: 1024)
# nze - number of non-zero entries (default: 128)

function sparseSignalGenerator(;n=4096,k=1,nze=128)
    A = spzeros(n,k)
    while nnz(A) < nze
        A[rand(1:n),rand(1:k)] = rand(Normal())
    end
    return A
end

# Function: arrayDecomposition
# Description: construct a 2nX1 vector using a nX1 vector following the following decomposition
    
#------------Inputs------------
# x0 - original sparse signal

function arrayDecomposition(x)
    
    l = Float64[]
    z = Float64[]
    
    for i in 1:length(x)
        push!(l, max(0,x[i]))
        push!(z, max(0,-x[i]))
    end

    return vcat(l,z)
end

function arrayComposition(x)
    u = Float64[]
    v = Float64[]

    n = Int(length(x)/2)

    for i in 1:n
        push!(u,x[i])
        push!(v,x[n+i])
    end

    return u-v
end

# Function: signalMonotoneNonlinearSystem

#----------Inputs----------
# x - an element of the n-dimensional euclidean space
# k - number of rows of the linear operator matrix (default: 1024)
# n - number of columns of the linear operator matrix (default: 4096)

function signalMonotoneNonlinearSystem(x,A)

    z = Float64[] # z = F(u) = min{u, Mu + r}
    
    # see [1] to understand the following lines of the code
    u = arrayDecomposition(x)
    b = A*x
    w =  0.1*norm(b,Inf)
    r = w*ones(length(x)) + vcat(-b,b)
    M=hcat(vcat(A,-A),vcat(A,-A))
    v = M*u + r

    for i in 1:length(u)
        push!(z,min(u[i],v[i])) # see [2]
    end

    return z, w

end

# adaptated from [3] compose right after the F is evaluated
function signalAlgorithm(x0, y, m, F, proj; ξ = 0.5, σ = 1.0e-4, ρ = 0.74, η = 0.5, θ = 0.5, ϵ = 1.0e-4)

    # inicializing variables
    k = 0;
    x = x0;
    newx = x0;
    d = 0;
    tk = ξ;
    zk = 0;
    t0 = time()
    perror = 0;

    meritvalue = 1

    while meritvalue > ϵ && k < 100

        # descent direction
        if k == 0
            d = -F(x0, y)[1]; d = arrayComposition(d)
        else
            # variables saving-calculations purposes only
            vec_one = F(newx, y)[1]
            vec_two = F(x,y)[1]
            vec_one = arrayComposition(vec_one)
            vec_two = arrayComposition(vec_two)
            norm_one = norm(vec_two)

            # continuing...
            γk = vec_one - vec_two
            λk = 1 + max(0,-dot(γk,tk*d)/norm(tk*d)^2)/norm_one
            yk = γk+λk*tk*norm_one*d
            sk = tk*d
            τa = norm(yk)^2/dot(sk,yk)  
            τb = dot(sk,yk)/norm(sk)^2
            τk = θ*τa + (1-θ)*τb
            βτk = dot(vec_one,yk)/dot(d,yk)-(τk+τa-τb)*(dot(vec_one,sk)/dot(d,yk))
            βk = max(βτk,η*(dot(vec_one,d)/norm(d)^2))
            d = -vec_one + βk*d

            x = newx
            tk = ξ
            
        end

        # tk determination
        if k == 0
            while dot(d+tk*d,d) < σ*tk*norm(d)^2
                tk = ρ*tk
            end
        else
            while dot(-vec_one+tk*d,d) < σ*tk*norm(d)^2
                tk = ρ*tk
            end
        end

        zk = newx + tk*d

        vec_three = F(zk,y)[1]
        norm_two = norm(vec_three)
        array_one = arrayComposition(vec_three)

        if meritvalue < ϵ   
            et = time() - t0
            perror = 1;
            k += 1
            return k, et, zk, array_one, meritvalue, perror
        else
            αk = dot(array_one,newx-zk)/norm_two^2
            newx = proj(x-αk*array_one)
        end

        w1 = F(x,y)[2]
        array_two,w2 = F(newx,y)
        array_two = arrayComposition(array_two)

        f = w2*norm(newx,1)+0.5*norm(y*newx-m,2)^2
        g = w1*norm(x,1)+0.5*norm(y*x-m,2)^2

        meritvalue = (f-g)/g

        k += 1

    end

    array_two = arrayComposition(F(newx,y)[1])
    
    et = time() - t0
    return k, et, newx, array_two, meritvalue, perror

end

# adaptated from [3] compose as later as possible
function signAlgorithm(x0, A, m, F, proj; ξ = 0.5, σ = 1.0e-4, ρ = 0.74, η = 0.5, θ = 0.5, ϵ = 1.0e-4)

    # inicializing variables
    k = 0;
    x = x0;
    x = arrayDecomposition(x)
    newx = x;
    d = 0;
    tk = ξ;
    zk = 0;
    t0 = time()
    perror = 0;

    meritvalue = 1

    while meritvalue > ϵ && k < 100

        # descent direction
        if k == 0
            d = -F(x0, A)[1]
        else
            # variables saving-calculations purposes only
            vec_one = F(newx, A)[1]
            vec_two = F(x,A)[1]
            norm_one = norm(vec_two)

            # continuing...
            γk = vec_one - vec_two
            λk = 1 + max(0,-dot(γk,tk*d)/norm(tk*d)^2)/norm_one
            yk = γk+λk*tk*norm_one*d
            sk = tk*d
            τa = norm(yk)^2/dot(sk,yk)  
            τb = dot(sk,yk)/norm(sk)^2
            τk = θ*τa + (1-θ)*τb
            βτk = dot(vec_one,yk)/dot(d,yk)-(τk+τa-τb)*(dot(vec_one,sk)/dot(d,yk))
            βk = max(βτk,η*(dot(vec_one,d)/norm(d)^2))
            d = -vec_one + βk*d

            x = newx
            x = arrayDecomposition(x)
            tk = ξ
            
        end

        # tk determination
        if k == 0
            while dot(d+tk*d,d) < σ*tk*norm(d)^2
                tk = ρ*tk
            end
        else
            while dot(-vec_one+tk*d,d) < σ*tk*norm(d)^2
                tk = ρ*tk
            end
        end

        zk = newx + tk*d
        newzk = arrayComposition(zk)

        vec_three = F(newzk,A)[1]
        norm_two = norm(vec_three)
        array_one = arrayComposition(vec_three)

        if meritvalue < ϵ   
            et = time() - t0
            perror = 1;
            k += 1
            return k, et, zk, array_one, meritvalue, perror
        else
            αk = dot(vec_three,newx-zk)/norm_two^2
            newx = proj(x-αk*vec_three)
        end

        x = arrayComposition(x)
        w1 = F(x,A)[2]
        newx = arrayComposition(newx)
        array_two,w2 = F(newx,y)
        array_two = arrayComposition(array_two)

        f = w2*norm(newx,1)+0.5*norm(y*newx-m,2)^2
        g = w1*norm(x,1)+0.5*norm(y*x-m,2)^2

        meritvalue = (f-g)/g

        k += 1

    end

    array_two = arrayComposition(F(newx,A)[1])
    
    et = time() - t0
    return k, et, newx, array_two, meritvalue, perror

end

# references

# [1] Figueiredo, M.A.; Nowak, R.D.; Wright, S.J. Gradient projection for sparse reconstruction,
# application to compressed sensing and other inverse problems. IEEE J. Sel. Top. Signal Process.
# 2007, 1, 586-597.

# [2] Xiao, Y.H.; Wang, Q. Y.; Hu, Q.J. Non-smooth equations based method for l1-norm problems
# with applications to compressed sensing. Nonlinear Anal. TMA 2011, 74, 3570-3577.

# [3] Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex
# constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438
