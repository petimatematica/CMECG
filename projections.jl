# Reference: Izmailov, A., & Solodiv, M. (2007). Otimização, volume 2: métodos computacionais. IMPA

id(x) = x # for unconstrained problems

projRplus(x)=max.(0.0,x) # non negative constraint

projbox(x; a=0.0, b=1.0) = min.(max.(a, x), b) # n dimensional box constraint
