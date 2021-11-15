using Plots, Interpolations, Optims

## The social planner problem we are going to solve here is
#  V(k, z) = max_{c, kp, h} log(Ak^α*h^(1-α)) - B*h + β*E[V(kp, zp)]
#       s.t. h ∈ [0, 1]

# Predefine the parameters of the model
const β = 0.98
const B = 0.5
const α = 2/3
const A = 1.0
const u(c, h) = log(c) - B*h^2.0

## First we solve a deterministic model (E[⋅], z, and zp can be dropped)
## Set solver parameters
const N_k  = 200         # Grid size of state vector
const N_h  = 80          # Grid size of decision vector for hours worked
kmax = A^(1/(1-α))       # Maximum sustainable capital stock, see TA note week 2
kmin = 0.0001
K = range(kmin, kmax, length=N_k)   # These are the states we are going
                                    #     to solve for
const H = range(0.0001, 1.0, length=N_h)  # Discrete set of possible labor decisions

## Run value-function iteration
V_guess = zeros(N_k) # Creates a vector of current guess what the value
                     # function looks like of size N_k
V_new   = copy(V_guess)    # Initialize these
optim_kp_1 = copy(V_guess)
optim_h_1  = copy(V_guess)

tol = 1.0e-2   # Greatest tolerance we will accept
err = 10.0*tol # Initialize an error greater than the tolerance
ii  = 1        # Initialize a counter
ii_max = 500   # Set a cap how many iterations we run (just to stop the solver
               # in case there is a coding error)
while ii < ii_max && err > tol
    println("i=$ii, err=$err")
    # For each k_j ∈ K
    for jj in 1:N_k
        k_j = K[jj]
        V_of_curr_k = zeros(N_k, N_h)
        # For current choice of k_j, try all combinations of kp and h and eval.
        #   them. Save that in a matrix
        for ll in 1:N_k
            kp = K[ll]
            for mm in 1:N_h
                h = H[mm]
                c = A*k_j^α*h^(1-α) - kp # Using c = f(k) - kp
                if c < 0.0001 # Consuming too little
                    V_of_curr_k[ll, mm] = -1.0e6
                else
                    V_of_curr_k[ll, mm] = u(c, h) + β*V_guess[ll]
                end
            end
        end
        (V_max, coordinate) = findmax(V_of_curr_k)
        ll = coordinate[1]; mm = coordinate[2]
        V_new[jj] = V_max
        optim_kp_1[jj] = K[ll]
        optim_h_1[ jj] = H[mm]
    end
    global err = maximum(maximum(abs.(V_new .- V_guess)))
    global ii += 1
    global V_guess = copy(V_new)
end

fig1 =
plot(K, V_guess)
display(fig1)

fig2 =
plot( K, optim_kp_1, label="k' given k")
plot!(K, optim_h_1,  label="h given k")
display(fig2)

## Second, we solve a stochastic model

#  Define stochastic component and expectation operator
const Z = [-0.2,  -0.1,  0.0,  0.2]
const P = [0.25, 0.25, 0.25, 0.25]
const N_z = length(Z)

function Exp(i_kp, V)
    a_sum = sum(V[i_kp, :].*P)
    return(a_sum)
end

## Set solver parameters
kmax = (exp(Z[end])*A)^(1/(1-α))    # Maximum sustainable capital stock, see TA note week 2
kmin = 0.0001
K = range(kmin, kmax, length=N_k)   # These are the states we are going
                                    #     to solve for

V_guess = zeros(N_k, N_z)
V_new   = copy(V_guess)  # Initialize these
optim_kp_1 = copy(V_guess)
optim_h_1  = copy(V_guess)

## Run value-function iteration
tol = 1.0e-2   # Greatest tolerance we will accept
err = 10.0*tol # Initialize an error greater than the tolerance
ii  = 1        # Initialize a counter
ii_max = 500   # Set a cap how many iterations we run (just to stop the solver
               # in case there is a coding error)
while ii < ii_max && err > tol
    println("i=$ii, err=$err")
    # For each k_j ∈ K and z_l ∈ Z
    for jj in 1:N_k
    for ll in 1:N_z
        k_j = K[jj]
        z_l = Z[ll]
        V_of_curr_k = zeros(N_k, N_h)
        # For current choice of (k_j, z_l), try all combinations of kp and h
        #   and eval. them. Save that in a matrix
        for mm in 1:N_k
            kp = K[mm]
            for nn in 1:N_h
                h = H[nn]
                c = exp(z_l)*A*k_j^α*h^(1-α) - kp # Using c = f(k) - kp
                if c < 0.0001 # Consuming too little
                    V_of_curr_k[mm, nn] = -1.0e6
                else
                    V_of_curr_k[mm, nn] = u(c, h) + β*Exp(mm, V_guess)
                end
            end
        end
        (V_max, coordinate) = findmax(V_of_curr_k)
        mm = coordinate[1]; nn = coordinate[2]
        V_new[jj, ll] = V_max
        optim_kp_1[jj, ll] = K[mm]
        optim_h_1[ jj, ll] = H[nn]
    end
    end
    global err = maximum(maximum(abs.(V_new .- V_guess)))
    global ii += 1
    global V_guess = copy(V_new)
end

fig3 =
plot( K[2:end], V_guess[2:end, 1], label="z_t = z_1")
plot!(K[2:end], V_guess[2:end, 2], label="z_t = z_2")
plot!(K[2:end], V_guess[2:end, 3], label="z_t = z_3")
plot!(K[2:end], V_guess[2:end, 4], label="z_t = z_4")
display(fig3)

fig4 =
plot( K, optim_h_1[:, 1])
plot!(K, optim_h_1[:, 2])
plot!(K, optim_h_1[:, 3])
plot!(K, optim_h_1[:, 4])
plot!(K, optim_kp_1[:, 1])
plot!(K, optim_kp_1[:, 2])
plot!(K, optim_kp_1[:, 3])
plot!(K, optim_kp_1[:, 4])
display(fig4)
