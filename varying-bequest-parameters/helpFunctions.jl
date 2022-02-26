function Expect(tp, CoHLI_tp, V)

        # CoHLI Cash-on-Hand Less Income
        EV = sum(V(CoHLI_tp + y_tp)*Pϵ[tp, i_ytp] for (i_ytp, y_tp) in enumerate(Y[tp, :]))
        return EV
end

function age2period(age; start_age=min_age)
  temp_ages = all_ages[all_ages .>= start_age]
  if minimum(temp_ages) > age
    return(missing)
  else
    period = findfirst(temp_ages .>= age)
  end
end

# From Github robertdkirkby / economics-with-julia
using Distributions
function tauchenmethod(mew, sigmasq, rho, znum, q)
  #Create states vector and transition matrix for the discrete markov process
  #   approximation of AR(1) process z'=rho*z+e, e~N(mew,sigmasq), by
  #   Tauchens method
  sigma=sqrt(sigmasq); #stddev of e
  zstar=mew/(1-rho); #expected value of z
  sigmaz=sigma/sqrt(1-rho^2); #stddev of z

  z=zstar*ones(znum,1) + range(-q*sigmaz, q*sigmaz, length=znum);
  omega=z[2]-z[1]; #Note that all the points are equidistant by construction.

  zi=z*ones(1,znum);
  zj=ones(znum,1)*z';

  P_part1=cdf.(Normal(), ((zj .+ omega/2 .- rho*zi) .- mew)./sigma);
  P_part2=cdf.(Normal(), ((zj .- omega/2 .- rho*zi) .- mew)./sigma);

  P=P_part1 .- P_part2;
  P[:,1]    = P_part1[:,1];
  P[:,znum] = 1 .- P_part2[:,znum];

  return z, P
end

function tick()
    t0 = time_ns()
    task_local_storage(:TIMERS, (t0, get(task_local_storage(), :TIMERS, ())))
    #@info " started timer at: $(now())"
end

function printTimeLeft(periods_left, periods_to_solve, time_in)
  base_T = Int(floor(log10(T)))+1
  str2print1 = sprintf1("Solving for period %$base_T.0i, ", periods_left)
  est_time_left = time_in/(periods_to_solve-periods_left)*periods_left
  if est_time_left > 60.0
    str2print2 = sprintf1("%4.2f min. remaining", est_time_left/60.0)
  else
    str2print2 = sprintf1("%4.1f sec. remaining", est_time_left)
  end
  printfmt(str2print1)
  printfmt(str2print2)
  flush(stdout)
  return(str2print1*str2print2)
end
