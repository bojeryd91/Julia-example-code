# Set parameters
const r  = 0.03
const Y_size = 9
const γ = 2
const β = 0.90
const u(c) = c > 0.0 ? c^(1.0-γ)/(1.0-γ) : -1.0e10

# Average income profile
annual_y_bars = [148200, 161589, 174355, 186496, 198014, 208908, 219178,
                 228825, 237847, 246246, 254021, 261172, 267699, 273603,
                 278882, 283538, 287570, 290978, 293762, 295923, 297460,
                 298373, 298662, 298327, 297368, 295786, 293580, 290750,
                 287296, 283218, 278516, 273191, 267242, 260669, 253472,
                 245652, 237207, 228139, 218447, 208131]./1.0e4
annual_y_bars = 57.07/4.0.*annual_y_bars./mean(annual_y_bars)
ages = 25:64.1
temp_interp = LinearInterpolation((ages), log.(annual_y_bars),
                                    extrapolation_bc=Line())
const min_age = 30
start_retirement = 65
const max_age = 85

work_ages = min_age:start_retirement-0.001
retirement_ages = start_retirement:max_age
all_ages  = cat(work_ages, retirement_ages, dims=1)
y_bars = exp.(temp_interp(work_ages))

## Income process Y
# Create work-life transitions   μ    σ^2  ρ    n           max # σ
(z_work, P_work) = tauchenmethod(-0.2, 0.15, 0.0, Y_size, 1.5)
P_work = repeat(transpose(P_work[1, :]), outer=length(work_ages))
Y_work = transpose(reduce(hcat, [y.*exp.(z_work)    for y in y_bars]))
# And for retiree life (no transitions)
#   https://www.pensionsmyndigheten.se/orange-bloggen/orange-bloggen/hur-mycket-av-din-slutlon-kan-du-rakna-med-att-fa-ut-i-pension
y_during_retirement = 0.7*sum(Y_work[end, :].*P_work[end, :])
Y_retiree = reshape([y_during_retirement for y in 1:Y_size
                                         for t in retirement_ages],
                                         length(retirement_ages), Y_size)
P_retiree = reshape([1.0/Y_size for i in 1:Y_size
                                   for t in retirement_ages],
                     length(retirement_ages), Y_size)

# Put it together
const Pϵ = cat(dims=1, P_work, P_retiree)
const Y  = cat(dims=1, Y_work, Y_retiree)

const T  = size(Y, 1)

## Initialize state grids D, L, X, and help grid C

# Lowest values for some choice variables
const c_min = 1.0e-3
const s_min = -10.0

# Grid X
const x_min = s_min < 0.0 ? s_min*(1.0+r) : s_min
const x_max = 500.0
x_grid1 = range(x_min, 0.9*x_min, length=20)
x_grid2 = range(x_min, x_max, length=200)
const X = sort(unique(cat(x_grid1, x_grid2, dims=1)))
const X_size = length(X)

# Grid S
const s_max  = x_max
const S      = copy(X)
const S_size = length(S)
