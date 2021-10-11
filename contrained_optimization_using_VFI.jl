using Optim, Plots, BenchmarkTools, Interpolations

## The social planner problem we are going to solve here is
#  V(k, z) = max_{c, kp, h} log(Ak^α*h^(1-α)) - B*h + β*E[V(kp, zp)]
#       s.t. h ∈ [0, 1]
## First we solve a deterministic model (E[⋅], z, and zp can be dropped)

## Predefine the parameters of the model
β = 0.98
B = 0.0
α = 2/3
A = 1.0

## Set solver parameters
N_k  = 100         # Grid size of state vector
N_c  = 500         # Grid size of decision vector for consumption
N_h  = 30          # Grid size of decision vector for hours worked
kmax = A^(1/(1-α)) # Maximum sustainable capital stock, see TA note week 2
kmin = 0.0001
K = range(kmin, kmax, length=N_k)   # These are the states we are going
                                    # to solve for
C = range(0.0001, kmax, length=N_c) # Discrete set of possible consumption dec.
H = range(0.0, 1.0, length=N_h)     # Discrete set of possible labor decisions

V_guess = [0.0 for k in K] # Creates a vector of current guess what the value
                           # function looks like of size N_k
V_prev  = copy(V_guess)    # Initialize these
optim_c = copy(V_guess)
optim_h = copy(V_guess)

u(c, h) = log(c) - B*h

## Run value-function iteration
tol = 1.0e-2   # Greatest tolerance we will accept
err = 10.0*tol # Initialize an error greater than the tolerance
ii  = 1        # Initialize a counter
ii_max = 5     # Set a cap how many iterations we run (just to stop the solver
               # in case there is a coding error)
while ii < 500 && err > tol
    println("i=$ii, err=$err")
    # At the beginning of each iteration, create an interpolation object
    V_interp = LinearInterpolation(K, V_guess)
    global V_prev = copy(V_guess)
    # For each k_j ∈ K
    for jj in 1:N_k
        k_j = K[jj]
        V_of_curr_k = zeros(N_c, N_h)
        # For current choice of k_j, try all combinations of c and h and eval.
        #   them. Save that in a matrix
        for ll in 1:N_c
            for mm in 1:N_h
                kp = A*k_j^α*H[mm]^(1-α) - C[ll] # Using kp = f(k) - c
                if kp < kmin # Saving too little, infeasible
                    V_of_curr_k[ll, mm] = -1.0e6
                else
                    V_of_curr_k[ll, mm] = u(C[ll], H[mm]) + β*V_interp(kp)
                end
            end
        end
        (_, coordinate) = findmax(V_of_curr_k)
        ll = coordinate[1]; mm = coordinate[2]
        V_guess[jj] = V_of_curr_k[ll, mm]
        optim_c[jj] = C[ll]
        optim_h[jj] = H[mm]
    end
    global err = maximum(maximum(abs.(V_guess .- V_prev)))
    global ii += 1
end

fig1 =
plot(K, V_guess)
display(fig1)

fig2 =
plot( K, optim_c, label="c given k")
plot!(K, optim_h, label="h given k")
