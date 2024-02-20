# Reference: Izmailov, A., & Solodiv, M. (2007). Otimização, volume 2: métodos computacionais. IMPA

function projRplus(x)
    
    v = Float64[]

    for i in 1:length(x)
        push!(v, max(0,x[i]))
    end

    return v

end

id(x) = x