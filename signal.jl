using LinearAlgebra, Distributions

# Function: sparseSignalGenerator
# Description: Generate a sparse matrix using the given values of number of rows, columns
# and non-zero entries, respectively.

#----------Inputs---------------
# n - number of rows (default: 4096)
# k - number of columns (default: 1024)
# nze - number of non-zero entries (default: 128)

function sparseSignalGenerator(;n=4096,k=1024,nze=128)
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

# Function: signalRecovery

#----------Inputs----------
# x0 - original sparse signal

function signalRecovery()
    x0 = sparseSignalGenerator() # original sparse signal
    y = rand(Normal(), size(x0)[2], size(x0)[1]) # linear operator
    m = y*x0

    x00 = y'*m  # guess

    gnoise = rand(Normal(0, 1.e-4), size(m)[1]) # Gaussian noise
    w = m + gnoise

    # WARNING: This code is not finished yet.

end