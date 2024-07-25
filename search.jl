# Reference: Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex
# constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438

function ding(x0, F, proj;xi=0.9,sigma=1.e-5,rho=0.74,eta=0.1,theta=0.1,maxiter=1000,tol=1.e-5)
    
    k=0
    t0 = time()
    Fevals = 0

    # step 1 (k=0)
    Fx0 = F(x0)
    Fevals +=1
    nFx0 = sqrt(Fx0'*Fx0)
    if nFx0 < tol
        et = time() - t0
        return x0, k, et, nFx0, Fevals, 0
    else
        dk = -Fx0
    end

    # step 2 (k=0)
    tk=xi
    ndk=sqrt(dk'*dk)
    e=sigma*ndk^2
    while -F(x0+tk*dk)'*dk < e*tk
        tk = rho*tk
        Fevals +=1
    end
    zk=x0+tk*dk

    # step 3 (k=0)
    Fzk=F(zk)
    Fevals +=1
    nFzk=sqrt(Fzk'*Fzk)
    if nFzk < tol
        et = time() - t0
        return x0, k, et, nFx0, Fevals, 1
    else
        alphak=(Fzk'*(x0-zk))/nFzk^2
        x1=proj(x0-alphak*Fzk)
        Fx1=F(x1)
        Fevals +=1
        nFx1=sqrt(Fx1'*Fx1)
    end

    # step 4 (k=0)
    k = k+1

    while k<maxiter
        
        # step 1 (k>0)
        if nFx1 < tol
            et = time() - t0
            return x1, k, et, nFx1, Fevals, 0
        else
            gammak=Fx1-Fx0
            sk = tk*dk
            nsk=sqrt(sk'*sk)
            lambdak=1.0+max(0.0,-(gammak'*sk)/nsk^2)/nFx0
            yk=gammak+lambdak*nFx0*sk
            c=sk'*yk
            d=dk'*yk
            tauA=(yk'*yk)/c
            tauB=c/nsk^2
            tauk=theta*tauA+(1.0-theta)*tauB
            betak=(Fx1'*yk)/d-(tauk+tauA-tauB)*(Fx1'*sk)/d
            betakplus=max(betak,eta*(Fx1'*dk)/ndk^2)
            dk = -Fx1+betakplus*dk
        end

        # step 2 (k>0)
        x0=x1
        Fx0=Fx1
        nFx0=nFx1

        tk=xi
        ndk=sqrt(dk'*dk)
        e=sigma*ndk^2
        while -F(x0+tk*dk)'*dk < e*tk
            tk = rho*tk
            Fevals +=1
        end
        zk=x0+tk*dk

        # step 3 (k>0)
        Fzk=F(zk)
        Fevals +=1
        nFzk=sqrt(Fzk'*Fzk)
        if nFzk < tol
            et = time() - t0
            return x0, k, et, nFx0, Fevals, 1
        else
            alphak=(Fzk'*(x0-zk))/nFzk^2
            x1=proj(x0-alphak*Fzk)
            Fx1=F(x1)
            Fevals +=1
            nFx1=sqrt(Fx1'*Fx1)
        end

        # step 4 (k=0)
        k = k+1
    end

    et = time() - t0
    return x1, k, et, nFx1, Fevals, 2

end
