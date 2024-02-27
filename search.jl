# Reference: Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex
# constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438

function algorithm(x0, F, proj; ξ = 1.0, σ = 1.0e-4, ρ = 0.74, η = 0.5, θ = 0.5, ϵ = 1.0e-5)

    # inicializing variables
    k = 0;
    x = x0;
    newx = x0;
    d = 0;
    tk = ξ;
    zk = 0;
    t0 = time()
    perror = 0;

    while norm(F(newx)) > ϵ


        fk = meritfunction(newx)
        fkm1 = meritfunction(x)

       
        test = abs(fk-fkm1) / abs(fkm1)
       
        println("iter = $k   tk= $tk   test = $test")

        if test < ϵ && k > 0
            return k, 0.0, newx, F(newx), norm(F(newx)), perror
        end

        # descent direction
        if k == 0
            d = -F(x0)
        else
            γk = F(newx) - F(x)
            λk = 1 + max(0,-dot(γk,tk*d)/norm(tk*d)^2)/norm(F(x))
            yk = γk+λk*tk*norm(F(x))*d
            sk = tk*d
            τa = norm(yk)^2/dot(sk,yk)  
            τb = dot(sk,yk)/norm(sk)^2
            τk = θ*τa + (1-θ)*τb
            βτk = dot(F(newx),yk)/dot(d,yk)-(τk+τa-τb)*(dot(F(newx),sk)/dot(d,yk))
            βk = max(βτk,η*(dot(F(newx),d)/norm(d)^2))
            d = -F(newx) + βk*d

            x = newx
            tk = ξ
            
        end
        # tk determination
        while dot(-F(newx+tk*d),d) < σ*tk*norm(d)^2
           # println("tk = $tk")
            tk = ρ*tk
        end

        zk = newx + tk*d

        if norm(F(zk)) < ϵ       #colocar código de erro
            et = time() - t0
            perror = 1;
            k += 1
            return k, et, zk, F(zk), norm(F(zk)), perror
        else
            αk = dot(F(zk),newx-zk)/norm(F(zk))^2
            newx = proj(x-αk*F(zk))
        end

        k += 1

    end
    
    et = time() - t0
    return k, et, newx, F(newx), norm(F(newx)), perror

end



function meritfunction(z)

    n = size(z,1)
    u = z[1:Int(n/2)] 
    v = z[Int(n/2)+1:n]
    x = u - v

    return 0.1 * norm(x,1) + 0.5 * norm(A*x - b)^2

end

