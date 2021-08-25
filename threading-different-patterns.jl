using DataFrames, Optim, Base.Threads, BenchmarkTools

N = Int(1e5)
input = rand(N)

function usingAppend(x)
        dfs = Vector{DataFrame}(undef, length(x))
        @threads for i in eachindex(dfs)
                df = DataFrame(x = x[i], y = 0.0)
                y = (x[i]).^2
                f(x) = (y .- x[1]).^2.0
                df[1, :y] = optimize(f, -1.0, 1.0).minimizer
                dfs[i] = df
        end
        sims = DataFrame()
        [append!(sims, df) for df in dfs]
        return(sims)
end
usingAppend([0.5, 0.2])

function usingReduce(x)
        dfs = Vector{DataFrame}(undef, length(x))
        @threads for i in eachindex(dfs)
                df = DataFrame(x = x[i], y = 0.0)
                y = (x[i]).^2
                f(x) = (y .- x[1]).^2.0
                df[1, :y] = optimize(f, -1.0, 1.0).minimizer
                dfs[i] = df
        end
        sims = reduce(vcat, dfs)
        return(sims)
end
usingReduce([0.5, 0.2])

println("Do both commands return the same?")
usingReduce([0.5, 0.2]) == usingAppend([0.5, 0.2]) ?
                println("Yes") : println("No")

@btime usingAppend(input)
@btime usingReduce(input)
