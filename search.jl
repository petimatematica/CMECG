# Reference: Yanyun Ding, Yunhai Xiao & Jianwei Li (2017) A class of conjugate gradient methods for convex
# constrained monotone equations, Optimization, 66:12, 2309-2328, DOI:10.1080/02331934.2017.1372438

# function ding(x0, F, proj;xi=0.9,sigma=1.e-5,rho=0.74,eta=0.5,theta=0.1,maxiter=1000,tol=1.e-5)
#-----------------------------------------------
# function ding implements the projected conjugate gradient method proposed in the reference above to find zeros of a vector field

# Input Parameters
# =================
# x0 (vector{Float64}) current estimation to be a zero of the vector field.
# F (function) objective function.
# proj (function) maps a point x into its projection into a feasible set.

# Optional Input Parameters
# =========================
# xi (Float64) the parameter ξ used in the Step 2 of the Algorithm 2.1 of the reference (default = 0.9).
# sigma (Float64) the parameter σ used in the Step 0 of the Algorithm 2.1 of the reference (default = 1.e-5).
# rho (Float64) the parameter ρ used in the Step 0 of the Algorithm 2.1 of the reference (default = 0.74).
# eta (Float64) the parameter η used in the Step 0 of the Algorithm 2.1 of the reference (default = 0.5).
# theta (Float64) the parameter θ used in the Step 0 of the Algorithm 2.1 of the reference (default = 0.1).
# maxiter (Int) set the maximum number of iterations (default = 1000).
# tol (Float64) tolerance given. (default = 1e-5).

# Output
# ======
# x (vector{Float64}) current estimation to be a zero of the vector field.
# k (Int) number of iterations.
# t (Float64) elapsed time to solve the problem
# nFx (Float64) contain the value of the euclidean norm of the current estimation.
# Fevals (Int) number of function evaluations.
#= ierror (Int) the value stored in this variable tells the following messages: 
0 - ||F(x0)||<tol, 1 - ||F(zk)||<tol, 2 - the maximum number of iterations has been exceeded.=#

function ding(x0, F, proj;xi=0.9,sigma=1.e-5,rho=0.74,eta=0.1,theta=0.1,maxiter=1000,tol=1.e-5)
    
    k=0 #number of iterations
    t0 = time() # initial time
    Fevals = 0 # number of function evaluations

    # step 1 (k=0)
    Fx0 = F(x0)
    Fevals +=1
    nFx0 = sqrt(Fx0'*Fx0) # avoiding the use of LinearAlgebra package
    if nFx0 < tol  #test if the guess is close enough to be regarded as a zero of the vector field
        et = time() - t0 # elapsed time
        return x0, k, et, nFx0, Fevals, 0
    else
        dk = -Fx0 
    end

    # step 2 (k=0)
    tk=xi
    ndk=sqrt(dk'*dk)
    e=sigma*ndk^2 # avoiding unecessary multiplications in the while loop
    while -F(x0+tk*dk)'*dk < e*tk # linesearch
        tk = rho*tk # backtracking
        Fevals +=1
    end
    zk=x0+tk*dk

    # step 3 (k=0)
    Fzk=F(zk)
    Fevals +=1
    nFzk=sqrt(Fzk'*Fzk)
    if nFzk < tol #test if the x0 is close enough to be regarded as a zero of the vector field
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
            c=sk'*yk # avoiding calculate the same inner product twice
            d=dk'*yk # avoiding calculate the same inner product twice
            tauA=(yk'*yk)/c
            tauB=c/nsk^2
            tauk=theta*tauA+(1.0-theta)*tauB
            betak=(Fx1'*yk)/d-(tauk+tauA-tauB)*(Fx1'*sk)/d
            betakplus=max(betak,eta*(Fx1'*dk)/ndk^2)
            dk = -Fx1+betakplus*dk
        end

        # step 2 (k>0)
        x0=x1 # since dk was defined the information about the previous iteration is unecessary
        Fx0=Fx1
        nFx0=nFx1

        tk=xi
        ndk=sqrt(dk'*dk)
        e=sigma*ndk^2 # avoiding unecessary multiplications in the while loop
        while -F(x0+tk*dk)'*dk < e*tk # linesearch
            tk = rho*tk # backtracking
            Fevals +=1
        end
        zk=x0+tk*dk

        # step 3 (k>0)
        Fzk=F(zk)
        Fevals +=1
        nFzk=sqrt(Fzk'*Fzk)
        if nFzk < tol #test if the x0 is close enough to be regarded as a zero of the vector field
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

    et = time() - t0 # elapsed time
    return x1, k, et, nFx1, Fevals, 2 # the solution was not found is less than maxiter iterations

end
