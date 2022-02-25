using Plots
# See https://en.wikipedia.org/wiki/Lagrange_polynomial
function LagrangeExtrap(x, y, x_new)
        l_j = [prod( j != m ? (x_new - x[m])/(x[j] - x[m]) : 1.0
                         for m in eachindex(x)) for j in eachindex(x)]
        return(sum(y.*l_j))
end

## Illustration
x = vcat(range(0.0, 1.0, length=100))
x_extp = x[end] .+ vcat(range(0.0, 0.5, length=10))
y = x.^2
plot(x, y, label="")
y_extp1 = LagrangeExtrap.((x[end-1:end],), (y[end-1:end],), x_extp)
y_extp2 = LagrangeExtrap.((x[end-2:end],), (y[end-2:end],), x_extp)
plot!(x_extp, y_extp1, label="")
plot!(x_extp, y_extp2, label="")

## Illustration of sinus example
f(x) = sin.(π*(x.-0.5))
sinX = f(x)
plot( x, sinX, label="", color=:blue)
plot!(x_extp, f.(x_extp), label="Actual f(x)", color=:blue,
        legend=:topleft)
plot!(x_extp, LagrangeExtrap.((x[end-2:end],), (sinX[end-2:end],), x_extp),
        label="2nd order approx.", linewidth=3,
        color=:red)
plot!(x_extp, LagrangeExtrap.((x[end-4:end],), (sinX[end-4:end],), x_extp),
        label="4th order approx.", linewidth=2,
        color=:green, linestyle=:dash)
