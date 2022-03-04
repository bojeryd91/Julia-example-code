using NLsolve, BenchmarkTools

function example1()
    f(X) = X.^3 .- A
    res = nlsolve(f, zeros(3))
    return(res.zero)
end

function example2(A_in)
    f(X) = X.^3 .- A_in
    res = nlsolve(f, zeros(3))
    return(res.zero)
end

function example3()
    f(X) = X.^3 .- A_const
    res = nlsolve(f, zeros(3))
    return(res.zero)
end

function simExample1(A_in)
    global A = A_in
    example1()
end

function simExample2(A_in)
    example2(A_in)
end

function simExample3()
    A = rand(3) # This line is only to add computation time
    example3()
end

##
println(""); println("----------------------------------------------")
A_same = rand(3)
const A_const = copy(A_same)
simExample1(A_same); simExample2(A_same); simExample3() # pre compiling
@btime simExample1(A_same)
@btime simExample2(A_same)
@btime simExample1(A_same)
@btime simExample2(A_same) # Passing as arguments seems to be quicker
@btime simExample3()
@btime simExample3()       # Seems that using const is as fast as passing
@btime simExample3()       # the variable as an argument
