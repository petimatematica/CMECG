# Reference:
# Surjanovic, S. & Bingham, D. (2013). Virtual Library of Simulation Experiments: Test Functions and Datasets. 
# Retrieved May 21, 2024, from http://www.sfu.ca/~ssurjano.

# Many Local Minima

# Ackley Funtion

function ackley(x; a=20, b=0.2, c=2*pi)
    d = size(x)[1] # calculates how many coordinates the vector x has
    s1 = sum(x[i]^2 for i in 1:d) # s1 means sum 1
    s2 = sum(cos(c*x[i]) for i in 1:d) # s2 means sum 2
    return -a*exp(-b*sqrt(s1/d))-exp(s2/d)+a+exp(1)
end

# Griewank Function
function griewank(x)
    d = size(x)[1] # calculates how many coordinates the vector x has
    s = sum(x[i]^2*0.00025 for i in 1:d) # s means sum
    p = prod(cos(x[i]/sqrt(i)) for i in 1:d) # p means product
    return s-p+1
end

# Levy Function
function levy(x)
    w = (x.-1)*0.25.+1
    d = size(x)[1] # calculates how many coordinates the vector x has
    s = sum((w[i]-1)^2*(1+10*sin(pi*w[i]+1)^2) for i in 1:(d-1))
    return sin(pi*w[1])+s+(w[d]-1)^2*(1+sin(2*pi*w[d])^2)
end

# Bowl-Shaped

# Perm Function 0, D, Beta

function perm(x; β=1)
    d = size(x)[1] # calculates how many coordinates the vector x has
    s = 0
    a = 0
    for i in 1:d
        for j in 1:d
            a+=(j+β)*(x[j]^i-1/j^i) # a means addend
        end
        s+=a^2
        a=0
    end
    return s
end

# Rotated Hyper-Ellipsoid Function

function rotated(x)
    d = size(x)[1] # calculates how many coordinates the vector x has
    s = 0
    a = 0
    for i in 1:d
        for j in 1:i
            a+=x[j]^2 # a means addend
        end
        s+=a
        a=0
    end
    return s
end

# Sphere Function

sphere(x)=sum(x[i]^2 for i in 1:size(x)[1])

# Sum of Different Powers Function

sumdiffpowers(x)=sum(abs(x[i])^(i+1) for i in 1:size(x)[1])

# Sum Squares Function

sumsquares(x)=sum(i*x[i]^2 for i in 1:size(x)[1])

# Trid Function

function trid(x)
    d = size(x)[1] # calculates how many coordinates the vector x has
    s1 = sum((x[i]-1)^2 for i in 1:d) # s1 means sum 1
    s2 = sum(x[i]*x[i-1] for i in 2:d) # s2 means sum 2
    return s1 - s2
end

# Plate-Shaped

# Zakharov Function

function zakharov(x)
    d = size(x)[1]
    s1 = sum(x[i]^2 for i in 1:d) # s1 means sum 1
    s2 = sum(0.5*i*x[i] for i in 1:d) # s2 means sum 2
    return s1 + s2^2 + s2^2
end

# Valley Shaped

# Dixon Price Function

dixon(x)=(x[1]-1)^2+sum(i*(2*x[i]^2-x[i-1])^2 for i in 2:size(x)[1])

# Steep Ridges/Drops

# Michalewicz Function

michalewicz(x; m=10)=-sum(sin(x[i])*sin((i*x[i]^2)/pi)^(2*m) for i in 1:size(x)[1])

# other

# Perm Function D, Beta

function permD(x; β=1)
    d = size(x)[1] # calculates how many coordinates the vector x has
    s = 0
    a = 0
    for i in 1:d
        for j in 1:d
            a+=(j^i+β)*((x[j]/j)^i-1) # a means addend
        end
        s+=a^2
        a=0
    end
    return s
end

# Styblinsky-Tang Function

styblitang(x) = 0.5*sum(x[i]^4-16x[i]^2+5x[i] for i in 1:size(x)[1])
