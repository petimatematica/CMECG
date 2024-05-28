include("search.jl")

using Random, Plots, BenchmarkProfiles

function myresults(F, proj;dim=50, nproblems=10, factor=10)

    etavalues=[.1; .3; .5]
    T = rand(nproblems,3)
    A = rand(nproblems) # time matrix
    B = rand(nproblems) # iterations matrix
    C = rand(nproblems) # function evaluations matrix

    for i in 1:size(etavalues)[1]

        rng = MersenneTwister(1234)

        for j in 1:nproblems
            k, et, fcounter, xk, Fxk, normfxk, perror = solver(factor*rand(rng, dim), F, proj, η=etavalues[i])
            T[j,1]=et
            T[j,2]=k
            T[j,3]=fcounter
            #println("$et, $j") 
        end

        A = hcat(A,T[1:end,1])
        B = hcat(B,T[1:end,2])
        C = hcat(C,T[1:end,3])

    end

    A = A[1:end, 2:end]
    B = B[1:end, 2:end]
    C = C[1:end, 2:end]
    p = performance_profile(PlotsBackend(), A, ["η=0.1", "η=0.3", "η=0.5"], xlabel = "CPU ratio time")
    q = performance_profile(PlotsBackend(), B, ["η=0.1", "η=0.3", "η=0.5"], xlabel = "CPU ratio iterations")
    r = performance_profile(PlotsBackend(), C, ["η=0.1", "η=0.3", "η=0.5"], xlabel = "CPU ratio funcion evaluations")

    return A, B, C, p, q, r

end