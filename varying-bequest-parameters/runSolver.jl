function runSolver(additional_end_wealth, bequest_util)

        ########################################################################
        ### Initialize some stuff
        ########################################################################

        # Create value function and decision matrices
        V_t = -1.0e10.*ones(T, X_size)
        c_t =   c_min.*ones(T, X_size)

        ########################################################################
        ### Run last period
        ########################################################################

        @threads for (i_x, x_state) in collect(enumerate(X))
                rhs_bc = x_state + additional_end_wealth
                c_t[T, i_x] = c_min
                if rhs_bc + Ψ_1 > c_min
                        to_minimize(c_in) =
                                -u(c_in) - bequest_util(rhs_bc - c_in)
                        result = optimize(to_minimize, c_min,
                                                rhs_bc + Ψ_1, GoldenSection())

                        V_t[T, i_x] = -result.minimum
                        c_t[T, i_x] =  result.minimizer
                end
        end

        ########################################################################
        ### Then, start in period T-1, and iterate backwards
        ########################################################################
        tick();

        num_periods_to_solve = T-1
        periods_left = num_periods_to_solve
        @inbounds for t in [T-t for t in 1:num_periods_to_solve]
                time_printed =
                        printTimeLeft(periods_left, num_periods_to_solve,
                                                                peektimer())
                periods_left -= 1
                # Prepare new interpolation object for period t+1
                V_tp = LinearInterpolation(X, V_t[t+1, :])

                #  Precomputation step
                W_t = -1.0e10.*ones(S_size)
                @threads for (i_s, s_choice) in collect(enumerate(S))
                        # Check if choice doesn't violate the borrowing
                        #       constraint.
                        if s_choice >= s_min
                                # Compute Cash-on-Hand Less Income
                                CoHLI = (1+ir(s_choice))*s_choice
                                W_t[i_s] = β*Expect(t+1, CoHLI, V_tp)
                        end
                end
                W_interp = LinearInterpolation(S, W_t)

                #   Solve for each state
                @threads for (i_x, x_state) in collect(enumerate(X))
                        # Check borrowing constraint isn't violated at minimum
                        #       consumption
                        if x_state - c_min > s_min
                                rhs(s) = u(x_state - s) + W_interp(s)
                                res = optimize(s -> -rhs(s),
                                                s_min, x_state - c_min,
                                                GoldenSection())
                                V_t[t, i_x] = -res.minimum
                                tilde_s     =  res.minimizer
                                c_t[t, i_x] =  x_state - tilde_s
                        end
                end

                print(" ")
                print(repeat("\b", 1+length(time_printed)))
        end

        print("Complete!"*repeat(" ", 40)*"\n")
        #show(to)
        @printf("\nTotal time: %3.1f min\n", peektimer()/60)
        flush(stdout)
        return (V_t, c_t)
end
