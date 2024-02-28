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

# adaptated from [3] compose right after the F is evaluated
function signalalgorithm(x0, F, proj; ξ = 1.0, σ = 1.0e-4, ρ = 0.74, η = 0.5, θ = 0.5, ϵ = 1.0e-5)

    # inicializing variables
    k = 0;                      # k is the number of iterations
    xkminusone = x0;            # xkminusone is the previous term of the sequence (xk)
    xk = xkminusone;            # xk is the current term of the sequence (xk)
    d = xk;              
    tk = ξ;                     # tk is the steplength parameter
    zk = xk;                 
    t0 = time();               
    meritvalue = 1;

    while true

        Fxkminusone = F(xkminusone)
        Fxk = F(xk)
        c = norm(Fxkminusone)

        println("iter = $k   meritvalue= $meritvalue")

        # descent direction
        if k == 0
            d = -F(x0)
        else
            sk = tk*d
            v = norm(sk)
            γk = Fxk - Fxkminusone
            λk = 1 + max(0,-dot(γk,sk)/v^2)/c
            yk = γk+λk*tk*c*d
            τa = norm(yk)^2/dot(sk,yk)  
            τb = dot(sk,yk)/v^2
            τk = θ*τa + (1-θ)*τb
            βτk = dot(Fxk,yk)/dot(d,yk)-(τk+τa-τb)*(dot(Fxk,sk)/dot(d,yk))
            βk = max(βτk,η*(dot(Fxk,d)/norm(d)^2))
            d = -Fxk + βk*d

            xkminusone = xk
            tk = ξ

        end

        # tk determination
        p = norm(d)
        while dot(-F(xk+tk*d),d) < σ*tk*p^2
            tk = ρ*tk
        end

        zk = xk + tk*d
        Fzk = F(zk)
        nFzk = norm(Fzk)

        αk = dot(Fzk,xk-zk)/nFzk^2
        xk = proj(xk-αk*Fzk)

        meritfxk = meritfunction(xk)
        meritfxkminusone = meritfunction(xkminusone)

        meritvalue = abs(meritfxk- meritfxkminusone)/abs(meritfxkminusone)

        if meritvalue < ϵ       #colocar código de erro
            et = time() - t0
            k += 1
            return k, et, xk, meritvalue
        end

        k += 1

    end

end

function meritfunction(z)

    n = size(z,1)
    u = z[1:Int(n/2)] 
    v = z[Int(n/2)+1:n]
    x = u - v

    return 0.1 * norm(x,1) + 0.5 * norm(A*x - b)^2

end

# references

# [1] Figueiredo, M.A.; Nowak, R.D.; Wright, S.J. Gradient projection for sparse reconstruction,
# application to compressed sensing and other inverse problems. IEEE J. Sel. Top. Signal Process.
# 2007, 1, 586-597.

# [2] Xiao, Y.H.; Wang, Q. Y.; Hu, Q.J. Non-smooth equations based method for l1-norm problems
# with applications to compressed sensing. Nonlinear Anal. TMA 2011, 74, 3570-3577.

# [3] Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex
# constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438
