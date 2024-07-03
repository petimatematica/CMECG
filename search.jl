# Reference: Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex
# constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438

using LinearAlgebra

#step 0
function solver(x0, F, proj; ξ = 1.0, σ = 1.0e-4, ρ = 0.74, η = 0.5, θ = 0.5, ϵ = 1.0e-5)

    # inicializing variables
    k = 0;                      # k is the number of iterations (end of step 0)
    xkminusone = x0;            # xkminusone is the previous term of the sequence (xk)
    xk = xkminusone;            # xk is the current term of the sequence (xk)
    d = xk;                     # d is the direction vector
    tk = ξ;                     # tk is the steplength parameter
    zk = xk;                    # zk is the point that is projected in the hyperplane (for more details see the reference cited above)
    t0 = time();                # t0 is the initial time
    perror = 0;                 # if perror == 0 then the solution comes from the first stop condition and if perror == 1 then the second stop condition has achieved
    fcounter = 0;               # stores the number of function evaluations
    
    # step 1
    while norm(F(xk)) > ϵ       # first stop condition

        #println("iter = $k   norm= $c")

        # descent direction
        if k == 0
            d = -F(x0)
            fcounter += 1
        else
            
            Fxkminusone = F(xkminusone)
            Fxk = F(xk)
            c = norm(Fxkminusone)
            fcounter += 2
            
            # y_(k-1) determinations
            sk = tk*d
            v = norm(sk)
            γk = Fxk - Fxkminusone
            λk = 1 + max(0,-dot(γk,sk)/v^2)/c
            yk = γk+λk*tk*c*d

            # βτk determination
            τa = norm(yk)^2/dot(sk,yk)  
            τb = dot(sk,yk)/v^2
            τk = θ*τa + (1-θ)*τb
            βτk = dot(Fxk,yk)/dot(d,yk)-(τk+τa-τb)*(dot(Fxk,sk)/dot(d,yk))

            # βk determination
            βk = max(βτk,η*(dot(Fxk,d)/norm(d)^2))

            # dk determination
            d = -Fxk + βk*d

            xkminusone = xk
            tk = ξ

        end

        # tk determination
        p = norm(d)

        # step 2
        while dot(-F(xk+tk*d),d) < σ*tk*p^2
            tk = ρ*tk
            fcounter += 1
        end

        zk = xk + tk*d
        Fzk = F(zk)
        nFzk = norm(Fzk)
        fcounter += 1

        # step 3
        if nFzk < ϵ       #colocar código de erro
            et = time() - t0
            perror = 1;
            k += 1
            return k, et, fcounter, zk, Fzk, nFzk, perror
        else
            αk = dot(Fzk,xk-zk)/nFzk^2
            xk = proj(xk-αk*Fzk)
        end

        # step 4
        k += 1

    end
    
    et = time() - t0
    Fxk = F(xk)
    fcounter += 1
    return k, et, fcounter, xk, Fxk, norm(Fxk), perror

end
