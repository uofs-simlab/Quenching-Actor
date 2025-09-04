#!/usr/bin/env julia

using Logging
global_logger(SimpleLogger(stderr, Logging.Warn))

using fitzhughNagumo
using NonlinearSolve
using BenchmarkTools

function run_bisection_experiment(xs, θ, usmin, usmax, states; method_name="Auto")
    println("=" ^ 50)
    println("Testing with parameters:")
    println("  gg = 1")
    println("  xs = $xs")
    println("  θ = $θ")
    println("  u_min = $usmin")
    println("  u_max = $usmax")
    println("  Method = $method_name")
    println("=" ^ 50)
    
    target = states.slow.ψ
    println("Target ψ = $target")
    
    # Check if bracket is active (has sign change)
    u1 = fitzhughNagumo.F((xs, θ, usmin); states=states)
    u2 = fitzhughNagumo.F((xs, θ, usmax); states=states)
    println("F(u_min) - target = $(u1 - target)")
    println("F(u_max) - target = $(u2 - target)")
    println("Sign change exists: $((u1 - target) * (u2 - target) < 0)")
    
    if (u1 - target) * (u2 - target) >= 0
        println("WARNING: No sign change in bracket - root finding may fail!")
    end
    
    # Define the problem
    prob = IntervalNonlinearProblem(
        (u, p) -> fitzhughNagumo.F((p[1], p[2], u); states=states) - p[3],
        (usmin, usmax),
        [xs, θ, target]
    )
    
    # Time the solution
    println("\nTiming the root finding...")
    result = @timed begin
        if method_name == "Bisection"
            println("Attempting to use Bisection method...")
            sol = solve(prob, Bisection(); abstol=fitzhughNagumo.atol, reltol=fitzhughNagumo.rtol)
            println("Solver used: $(typeof(sol.alg))")
        else
            sol = solve(prob; abstol=fitzhughNagumo.atol, reltol=fitzhughNagumo.rtol)
            println("Solver used: $(typeof(sol.alg))")
        end
    end
    
    sol = result.value
    runtime = result.time
    
    println("Runtime: $(runtime) seconds")
    println("Success: $(NonlinearSolve.SciMLBase.successful_retcode(sol))")
    
    if NonlinearSolve.SciMLBase.successful_retcode(sol)
        println("Root found: u = $(sol.u)")
        function_at_root = fitzhughNagumo.F((xs, θ, sol.u); states=states)
        println("Function value at root: $(function_at_root)")
        println("Residual: $(abs(function_at_root - target))")
        
        # Check if solution is at boundary
        if isapprox(sol.u, usmin, atol=1e-10) || isapprox(sol.u, usmax, atol=1e-10)
            println("WARNING: Solution is at boundary!")
        end
        
        # Additional debugging - check function values around the solution
        println("\nDebugging function behavior around solution:")
        test_points = [sol.u - 1e-6, sol.u, sol.u + 1e-6]
        for (i, u_test) in enumerate(test_points)
            if u_test >= usmin && u_test <= usmax
                f_val = fitzhughNagumo.F((xs, θ, u_test); states=states)
                println("  u = $(u_test): F(u) = $(f_val), F(u) - target = $(f_val - target)")
            end
        end
        
        return sol.u, runtime, true
    else
        println("Root finding FAILED")
        println("Return code: $(sol.retcode)")
        return 1.0, runtime, false  # Return fail value like the original code
    end
end

function main()
    # Fixed parameters for the experiment
    xs = 34.2122
    θ = 5.0
    gg = 1
    
    # Set up file paths (same as your actor code)
    base = "/globalhome/tus210/HPC/quenchin_actor/"
    waved = base * "waves/index_11"
    crit = base * "waves/index_10"
    Ufile = waved * "/" * string(gg) * "/U"
    Pfile = waved * "/" * string(gg) * "/p"
    ufile = crit * "/" * string(gg) * "/U"
    pfile = crit * "/" * string(gg) * "/p"
    
    println("Loading states from files...")
    println("Ufile: $Ufile")
    println("Pfile: $Pfile")
    println("ufile: $ufile")
    println("pfile: $pfile")
    
    # Load states
    states = fitzhughNagumo.formStatesFromFiles(Ufile, Pfile, ufile, pfile)
    println("States loaded successfully!")
    println("Target ψ (slow.ψ): $(states.slow.ψ)")
    
    # Test case 1: Wide range [-1000, 0] with default solver
    println("\n" * "=" ^ 60)
    println("EXPERIMENT 1: Wide range [-1000, 0] with Default Solver")
    println("=" ^ 60)
    result1, time1, success1 = run_bisection_experiment(xs, θ, -1000.0, 0.0, states; method_name="Bisection")
    
    # Test case 2: Narrow range [-5.5, 0] with default solver
    println("\n" * "=" ^ 60)
    println("EXPERIMENT 2: Narrow range [-0.976562, 0] with Default Solver")
    println("=" ^ 60)
    result2, time2, success2 = run_bisection_experiment(xs, θ, -0.9765625, 0.0, states; method_name="Bisection")
    
    # Summary
    println("\n" * "=" ^ 60)
    println("EXPERIMENT SUMMARY")
    println("=" ^ 60)
    println("Parameters: gg=$gg, xs=$xs, θ=$θ")
    println()
    println("Test 1 ([-1000, 0] Default):")
    println("  Runtime: $(time1) seconds")
    println("  Success: $success1")
    println("  Result: $result1")
    println()
    println("Test 2 ([-5.5, 0] Default):")
    println("  Runtime: $(time2) seconds") 
    println("  Success: $success2")
    println("  Result: $result2")
    println()
    
    if success1 && success2
        speedup = time1 / time2
        println("Speedup factor (wide/narrow): $(speedup)x")
        println("Time difference: $(time1 - time2) seconds")
        println("Result difference: $(abs(result1 - result2))")
    end
    
    # Save results to file
    open("bisection_experiment_results.txt", "w") do io
        println(io, "Default Solver Root Finding Experiment Results")
        println(io, "Generated: $(Dates.now())")
        println(io, "Parameters: gg=$gg, xs=$xs, θ=$θ")
        println(io, "")
        println(io, "Test 1 - Wide range [-1000, 0] with Default Solver:")
        println(io, "  Runtime: $(time1) seconds")
        println(io, "  Success: $success1")
        println(io, "  Result: $result1")
        println(io, "")
        println(io, "Test 2 - Narrow range [-5.5, 0] with Default Solver:")
        println(io, "  Runtime: $(time2) seconds")
        println(io, "  Success: $success2") 
        println(io, "  Result: $result2")
        println(io, "")
        
        if success1 && success2
            speedup = time1 / time2
            println(io, "Speedup factor (wide/narrow): $(speedup)x")
            println(io, "Time difference: $(time1 - time2) seconds")
            println(io, "Result difference: $(abs(result1 - result2))")
        end
    end
    
    println("Results saved to: bisection_experiment_results.txt")
end

# Add Dates for timestamp
using Dates

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
