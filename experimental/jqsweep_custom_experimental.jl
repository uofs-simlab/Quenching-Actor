#!/usr/bin/env julia

using Logging
global_logger(SimpleLogger(stderr, Logging.Warn))

# Include the experimental fitzhughNagumo module with early termination support
include("fitzhughNagumo_experimental/src/fitzhughNagumo.jl")
using .fitzhughNagumo

"""
Custom fast bisection for quenching threshold detection
Supports two modes:
1. No early termination (original behavior)
2. Early termination using Sundials events (recommended for performance)
"""
function custom_G_with_early_termination(q; f=fitzhughNagumo.F, states, lims, fail=1.0, 
                                         max_iter=50, tolerance=1e-5, precomputed=nothing, 
                                         early_termination=true)
    xs, θ = q
    target = states.slow.ψ
    a_min, a_max = lims
    
    # Select the appropriate function based on early termination setting
    if early_termination
        f_func = fitzhughNagumo.F_with_sundials_events
    else
        f_func = fitzhughNagumo.F
    end
    
    # Initial function evaluations - use precomputed values if available
    if precomputed !== nothing && haskey(precomputed, a_min)
        f_min = precomputed[a_min] - target
    else
        f_min = f_func((xs, θ, a_min); states=states, early_termination=early_termination) - target
    end
    
    if precomputed !== nothing && haskey(precomputed, a_max)
        f_max = precomputed[a_max] - target
    else
        f_max = f_func((xs, θ, a_max); states=states, early_termination=early_termination) - target
    end
    
    # Check for sign change
    if f_min * f_max > 0
        return fail
    end
    
    # Ensure proper ordering (f_min < 0, f_max > 0)
    if f_min > 0
        a_min, a_max = a_max, a_min
        f_min, f_max = f_max, f_min
    end
    
    # Fast bisection loop with early termination support
    for iteration in 1:max_iter
        a_mid = (a_min + a_max) / 2
        f_mid = f_func((xs, θ, a_mid); states=states, early_termination=early_termination) - target
        
        # Check convergence
        if abs(f_mid) < tolerance || abs(a_max - a_min) < tolerance
            # Check if solution is at boundary (treat as failure)
            if abs(a_mid - lims[1]) < 1e-12 || abs(a_mid - lims[2]) < 1e-12
                return fail
            end
            return a_mid
        end
        
        # Update bracket
        if f_mid < 0
            a_min = a_mid
            f_min = f_mid
        else
            a_max = a_mid
            f_max = f_mid
        end
    end
    
    # Final check for boundary solution
    a_final = (a_min + a_max) / 2
    if abs(a_final - lims[1]) < 1e-12 || abs(a_final - lims[2]) < 1e-12
        return fail
    end
    
    return a_final
end

function main()
    if length(ARGS) < 10
        error("Usage: jqsweep_custom_experimental.jl Ufile Pfile ufile pfile n usmax gg xs theta early_termination")
    end

    # 1) parse arguments
    Ufile, Pfile, ufile, pfile = ARGS[1:4]
    n     = parse(Int64, ARGS[5])      
    usmax = parse(Float64, ARGS[6])
    gg    = parse(Float64, ARGS[7])
    xs    = parse(Float64, ARGS[8])
    θ     = parse(Float64, ARGS[9])
    early_termination = parse(Bool, ARGS[10])  

    # println("Parameters: xs=$xs, θ=$θ, early_termination=$early_termination, strategy=$termination_strategy")

    # 2) compute usmin from n: usmin = -1000/(2^n)
    usmin = -1000.0 / (2.0^n)

    # 3) load data into a States object
    states = fitzhughNagumo.formStatesFromFiles(Ufile, Pfile, ufile, pfile)

    # 4) build brackets (strategy from jqsweep_update.jl) - fixed to avoid degenerate brackets
    if n == 0
        # No bracket optimization - use single wide bracket like detailed_bisection_comparison.jl
        brackets = [(-1000.0, 0.0)]
    else
        # Bracket optimization enabled - use multiple brackets
        brackets = unique([
            (usmin, usmax),      # Optimized bracket
            (-1000.0, usmin)     # Fallback bracket
        ])
    end

    # Pre-check for sign change in widest bracket
    # Select function based on early termination setting
    if early_termination
        f_func = fitzhughNagumo.F_with_sundials_events
    else
        f_func = fitzhughNagumo.F
    end
    
    big_lo, big_hi = -1000.0, 0.0
    start_time = time()
    u_big_lo = f_func((xs, θ, big_lo); states=states, early_termination=early_termination) 
    boundary_lo_time = time() - start_time
    
    start_time = time()
    u_big_hi = f_func((xs, θ, big_hi); states=states, early_termination=early_termination)
    boundary_hi_time = time() - start_time
    
    
    target = states.slow.ψ
    if (u_big_lo - target) * (u_big_hi - target) >= 0
        # println("No sign change in bracket - returning failure")
        println("1.0")
        return
    end        
    
    # println("Sign change detected: ψ_lo=$(round(u_big_lo, digits=4)), ψ_hi=$(round(u_big_hi, digits=4)), target=$(round(target, digits=4))")        

    # Store the computed values to avoid recomputation in custom_G
    precomputed_vals = Dict(
        -1000.0 => u_big_lo,
        0.0 => u_big_hi
    )

    # 5) sweep until we get a non-failure A_try using CUSTOM bisection with early termination
    A = 1.0                # default = fail
    found_lo, found_hi = brackets[1]
  
    bisection_start_time = time()
    
    for (i, (u_lo, u_hi)) in enumerate(brackets)
        # println("Trying bracket $i: [$u_lo, $u_hi]")
        bracket_start_time = time()
        
        A_try = custom_G_with_early_termination((xs, θ);
                 states = states,
                 lims    = (u_lo, u_hi),
                 fail    = 1.0,
                 precomputed = precomputed_vals,
                 early_termination = early_termination)
        
        bracket_time = time() - bracket_start_time
        
        if A_try != 1.0
            A = A_try
            found_lo, found_hi = u_lo, u_hi
            # println("SUCCESS in bracket $i after $(round(bracket_time, digits=2))s: A* = $(round(A, digits=6))")
            break
        else
            # println("FAILED bracket $i after $(round(bracket_time, digits=2))s")
        end
    end
    
    total_bisection_time = time() - bisection_start_time
    # println("Total bisection time: $(round(total_bisection_time, digits=2))s")
    
    # if A != 1.0
    #     println("FINAL RESULT: A* = $(round(A, digits=8))")
    # else
    #     println("BISECTION FAILED - no valid solution found")
    # end

    # 6) emit A and the bracket that worked
    println("$A")
end

main()
