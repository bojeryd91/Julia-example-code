using Plots, Interpolations, Printf, TimerOutputs, Base.Threads, Optim,
        TickTock, Formatting, Plots.PlotMeasures

include("helpFunctions.jl")
include("initializeParametersAndGrids.jl")
include("runSolver.jl")

##
Ψ_0 = 50.0
Ψ_1 = 0.0

bequest_util(b) = b > 0.0 ? Ψ_0*(b+Ψ_1)^(1-γ)/(1-γ) : -1.0e10

add_wealth = 50.0
(V, C) = runSolver(add_wealth, bequest_util)

plot(X, C[end, :])
plot!([-Ψ_1], [0.0], seriestype=:scatter)

plot(X, X .- C[end-7, :])

##
C_itp = [LinearInterpolation(X, C[t, :], extrapolation_bc=Line()) for t in 1:T]

function df(t, x_t)
        beq = 0.0
        c_t = C_itp[t](x_t)
        s_t = x_t - c_t
        if t == T
                beq = x_t + add_wealth - c_t
                s_t = 0.0
        end
        return(c_t, s_t, beq)
end

##
x0 = 0.0
sim = zeros(T, 5)
sim[1, 2] = x0
for t in 1:T
        sim[t, 1] = t
        (c, s, beq) = df(t, sim[t, 2])
        sim[t, 3] = c
        sim[t, 4] = s
        if t < T # Compute next period's cash-on-hand
                sim[t+1, 2] = s*(1.0+r) + Y[t+1, 5]
        else     # Compute bequest
                sim[T, 5] = beq
        end
end
plot( sim[:, 1], sim[:, 2], color=:red,   label="")
plot!(sim[:, 1], sim[:, 3], color=:green, label="")
subplot = twinx()
plot!(subplot, [sim[end, 1]], [sim[end, 5]], color=:blue,
                label="", seriestype=:scatter, xaxis=false)
ylims!(subplot, 0.0, 1.1*sim[end, 5])
plot!(right_margin=10mm)
