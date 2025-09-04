#!/usr/bin/env julia

using Logging
global_logger(SimpleLogger(stderr, Logging.Warn))

# Use the fitzhughNagumo package (no include needed, it's in Project.toml)
using fitzhughNagumo

"""
Custom fast bisection for quenching threshold detection - optimized for actor workflows
"""
function custom_G(q; f=fitzhughNagumo.F, states, lims, fail=1.0, 
                  max_iter=50, tolerance=1e-5, precomputed=nothing)
    xs, θ = q
    target = states.slow.ψ
    a_min, a_max = lims
    
    # Initial function evaluations - use precomputed values if available
    if precomputed !== nothing && haskey(precomputed, a_min)
        f_min = precomputed[a_min] - target
    else
        f_min = f((xs, θ, a_min); states=states) - target
    end
    
    if precomputed !== nothing && haskey(precomputed, a_max)
        f_max = precomputed[a_max] - target
    else
        f_max = f((xs, θ, a_max); states=states) - target
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
    
    # Fast bisection loop
    for iteration in 1:max_iter
        a_mid = (a_min + a_max) / 2
        f_mid = f((xs, θ, a_mid); states=states) - target
        
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
    if length(ARGS) < 9
        error("Usage: jqsweep_custom.jl Ufile Pfile ufile pfile n usmax gg xs theta")
    end

    # 1) parse arguments
    Ufile, Pfile, ufile, pfile = ARGS[1:4]
    n     = parse(Int64, ARGS[5])      # Power for bracket computation
    usmax = parse(Float64, ARGS[6])
    gg    = parse(Float64, ARGS[7])
    xs    = parse(Float64, ARGS[8])
    θ     = parse(Float64, ARGS[9])

    # 2) compute usmin from n: usmin = -1000/(2^n)
    usmin = -1000.0 / (2.0^n)

    # 3) load states
    states = fitzhughNagumo.formStatesFromFiles(Ufile, Pfile, ufile, pfile)

    # 3) build bracket list: dynamic bracket first, then fallbacks
    # Optimize: we know umax is always 0.0, reuse usmin when it's -1000
    brackets = unique([
        (usmin, usmax),
        (-1000.0, usmin)
    ])

    # Pre-check for sign change in widest bracket
    # Optimize: if usmin is already -1000, reuse the computation
    big_lo, big_hi = -1000.0, 0.0
    u_big_lo= fitzhughNagumo.F((xs, θ, big_lo); states=states) 
    u_big_hi = fitzhughNagumo.F((xs, θ, big_hi); states=states)
    target = states.slow.ψ
    if (u_big_lo - target) * (u_big_hi - target) >= 0
        println("1.0")
        return
    end        

    # Store the computed values to avoid recomputation in custom_G
    precomputed_vals = Dict(
        -1000.0 => u_big_lo,
        0.0 => u_big_hi
    )

    # 4) sweep until we get a non-failure A_try using CUSTOM bisection
    A = 1.0                # default = fail
    found_lo, found_hi = brackets[1]
    for (u_lo, u_hi) in brackets
        # Use CUSTOM G function instead of the built-in one 
        # Pass precomputed values to avoid redundant evaluations
        A_try = custom_G((xs, θ);
                 states = states,
                 lims    = (u_lo, u_hi),
                 fail    = 1.0,
                 precomputed = precomputed_vals)
        if A_try != 1.0
            A = A_try
            found_lo, found_hi = u_lo, u_hi
            break
        end
    end

    # 5) emit A and the bracket that worked
    println("$A")
end

main()
