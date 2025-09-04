module fitzhughNagumo

using DelimitedFiles, StaticArrays, SparseArrays
using DiffEqOperators, DifferentialEquations, Sundials, LinearAlgebra
using DiffEqCallbacks: PeriodicCallback, FunctionCallingCallback, ContinuousCallback, DiscreteCallback, CallbackSet
using Dierckx
using NonlinearSolve
using Plots

export formStatesFromFiles, F, G, save_wave_plot, F_with_sundials_events
export h, atol, rtol

const Lglob = 2700.0
const N     = 1 + 2^13
const h     = Lglob / (N - 1)
const ξ     = range(-Lglob / 2, Lglob / 2, length = N)
const δ     = CenteredDifference(2, 6, h, N)
const B     = PeriodicBC(Float64)
const Δ     = sparse((δ * B))[1]
const atol  = 1e-8
const rtol  = 1e-6

function fhn!(dx, x, p, t)
  dx[1] = x[1] * (1 - x[1]) * (x[1] - p[2]) - x[2]
  dx[2] = p[3] * (p[1] * x[1] - x[2])
  return nothing
end

function fhn(x::SVector{2, T}, p, t) where T
  dx = MVector{2, T}(undef)
  fhn!(dx, x, p, t)
  return SVector(dx)
end

function fhn(x::MVector{2, T}, p, t) where T
  dx = MVector{2, T}(undef)
  fhn!(dx, x, p, t)
  return SVector(dx)
end


function fhn(x::AbstractVector{T}, p, t) where T
  dx = MVector{2, T}(undef)
  fhn!(dx, x, p, t)
  return SVector(dx)
end

function fhn_pde!(dx, x, p, t)
  @views du = dx[1, :]
  @inbounds for i in eachindex(du)
    fhn!(view(dx, :, i), view(x, :, i), p, t)
  end
  mul!(du, Δ, x[1, :], 1.0, 1.0)
  return nothing
end

function Ψ(u::AbstractVector{T}; ū=[zero(T), zero(T)]) where T<:Real
  spl = Spline1D(ξ, abs.(u .- ū[1]); k=5)
  return integrate(spl, minimum(ξ), maximum(ξ))
end

function Ψ(u::AbstractMatrix{T}; ū=[zero(T), zero(T)]) where T<:Real
  return sum(Ψ(row; ū=ū) for row in eachrow(u))
end

struct wave{T}
  u::Matrix{T}
  ψ::T
  p::Vector{T}
end

function rest(p::Vector{T}; u0 = SA{T}[0.0, 0.0]) where T
  fhn_func = (u, p, t) -> fhn(u, p, zero(T))
  problem = NonlinearProblem((u, p) -> fhn_func(u, p, 0), u0, p)
  sol     = solve(problem; abstol=zero(T), reltol=eps(T))
  return SciMLBase.successful_retcode(sol) ? sol.u : [NaN, NaN]
end

function save_wave_plot(U::Matrix, ξ::AbstractVector; filename="wave.png", title_str="Wave Profile",
                        label1="u₁", label2="u₂")
    p = plot(ξ, U[1, :], label=label1, lw=2, xlabel="ξ", ylabel="value",
             legend=:topright, title=title_str)
    plot!(p, ξ, U[2, :], label=label2, lw=2, linestyle=:dash)
    savefig(p, filename)
    # println("Saved plot to $filename")
end


function normalizeAutoSolution(M, L, ξ)
  x = M[:, 1]; u1 = M[:, 2]; u2 = M[:, 3]
  x0 = x[argmax(u1)]
  xtile = vcat(x, x .+ L)
  utile = vcat(hcat(u1, u2), hcat(u1, u2))
  xtile .-= x0
  s1 = Spline1D(xtile, utile[:, 1]; k = 5)
  s2 = Spline1D(xtile, utile[:, 2]; k = 5)
  U = zeros(eltype(x), 2, length(ξ))
  U[1, :] .= s1.(ξ)
  U[2, :] .= s2.(ξ)
  return U
end

function load_and_center(Ufile::String, Pfile::String)
  M   = readdlm(Ufile)
  pf  = vec(readdlm(Pfile))
  U   = normalizeAutoSolution(M, pf[5], ξ)
  p   = pf[1:3]
  ū   = rest(p)
  ψ   = Ψ(U[1, :]; ū=ū)
  return U, p, ψ, ū
end

function formStatesFromFiles(fastU::String, fastP::String, slowU::String, slowP::String)
  resolve(f) = isfile(f) ? f : joinpath(@__DIR__, f)
  fastU, fastP = resolve(fastU), resolve(fastP)
  slowU, slowP = resolve(slowU), resolve(slowP)
  Ufast, pfast, ψfast, ū  = load_and_center(fastU, fastP)
  Uslow, pslow, ψslow, _ = load_and_center(slowU, slowP)
  @info "Loaded Ψ_fast=ψfast   Ψ_slow=ψslo"
  pfull_fast = vec(readdlm(fastP))
  tspan = (0.0, Lglob / (2 * pfull_fast[4]))
  prob  = ODEProblem(fhn_pde!, deepcopy(Ufast), tspan, pfast)
  return (
    ū    = ū,
    fast = wave{Float64}(Ufast, ψfast, pfast),
    slow = wave{Float64}(Uslow, ψslow, pslow),
    prob = prob
  )
end


function X(q; H = sign)
  xs, θ, A = q
  out = zeros(typeof(xs), 2, N)
  @inbounds out[1, :] .= (A / 4) .* (1 .+ H.(ξ .- θ .+ xs / 2)) .* (1 .- H.(ξ .- θ .- xs / 2))
  return out
end

function F(q; states, early_termination=false)
  prob2 = remake(states.prob, u0 = states.prob.u0 .+ X(q))
  
  # println("F: Starting simulation with tspan = $(prob2.tspan)")
  # start_time = time()
  
  sol = solve(prob2, CVODE_BDF(linear_solver = :GMRES);
              abstol = atol, reltol = rtol,
              save_end = true, save_everystep = false)
  
  # end_time = time()
  # total_time = end_time - start_time
  
  final_ψ = Ψ(sol[1, :, end]; ū = states.ū)
  
  # println("F: Simulation completed in $(round(total_time, digits=3))s")
  # println("F: Final time reached: $(sol.t[end]) (target was $(prob2.tspan[2]))")
  # println("F: Final ψ = $(round(final_ψ, sigdigits=6))")
  # println("F: Solution status: $(sol.retcode)")
  
  return final_ψ
end

function F_with_sundials_events(q; states, early_termination=true, 
                               rest_threshold=0.01,     # threshold for max(u-0) for quenching
                               recovery_psi_tolerance=0.5,   # ψ must be within this distance of fast.ψ
                               recovery_deriv_threshold=1e-4,  # |ψ'| threshold for recovery stability
                               recovery_accel_threshold=1e-5,  # |ψ''| threshold for recovery stability
                               min_time_fraction=0.2,  # minimum simulation fraction before early exit
                               history_length=5)        # number of ψ values to keep for derivative computation
  """
  F function using Sundials native event detection for optimal performance.
  
  CORRECTED PHYSICS-BASED EVENT DETECTION:
  - Quenching (A < A*): Event triggered when max(|u-0|) < rest_threshold (simple check)
  - Recovery (A > A*): Event triggered when ψ ≈ fast.ψ AND |ψ'(t)| < tol AND |ψ''(t)| < tol
    The recovered wave stabilizes AT the fast wave integral, not above it!
  - At threshold (A ≈ A*): No events triggered, runs to completion → ψ ≈ ψ_slow
  """
  
  if !early_termination
    return F(q; states=states)
  end
  
  
  prob2 = remake(states.prob, u0 = states.prob.u0 .+ X(q))
  t0, t1 = prob2.tspan
  min_time = t0 + min_time_fraction * (t1 - t0)
  
  recovery_psi_min = states.fast.ψ - recovery_psi_tolerance
  recovery_psi_max = states.fast.ψ + recovery_psi_tolerance
  
  ψ_history = Float64[]
  t_history = Float64[]
  
  quench_condition(u, t, integrator) = begin
    # Only check after minimum time
    if t < min_time
      return 1.0  # Positive = no event
    end
    
    
    u_component = @view u[1, :]
    
    # QUENCHING DETECTION: Simple and fast
    max_u_deviation = maximum(abs.(u_component))
    if max_u_deviation < rest_threshold
      return -1.0  # Negative = trigger event (quenching detected)
    end
    
    return 1.0  # Positive = continue integration
  end
  
  # Quenching event action
  quench_affect!(integrator) = begin
    u_component = @view integrator.u[1, :]
    max_u = maximum(abs.(u_component))
    
    # println("SUNDIALS QUENCH EVENT @ t=$(round(integrator.t, sigdigits=4)): QUENCHING (A < A*)")
    # println("max(|u|) = $(round(max_u, sigdigits=4)) < rest_threshold($(rest_threshold))")
    # println("Progress: $(round(100 * (integrator.t - t0)/(t1 - t0), digits=1))% of simulation time")
    
    DiffEqBase.terminate!(integrator)
  end
  
  # Recovery detection using periodic callback (needs derivative computation)
  recovery_affect!(integrator) = begin
    # Only check after minimum time
    if integrator.t < min_time
      return
    end
    
    # Extract u component and compute ψ
    u_component = @view integrator.u[1, :]
    current_psi = Ψ(u_component; ū = states.ū)
    
    # Store history
    push!(ψ_history, current_psi)
    push!(t_history, integrator.t)
    
    # Keep only recent history
    if length(ψ_history) > history_length
      popfirst!(ψ_history)
      popfirst!(t_history)
    end
    
    # Need at least 3 points for second derivative
    if length(ψ_history) >= 3
      n = length(ψ_history)
      
      # Check if ψ is close to stable wave integral (fast.ψ)
      is_near_stable = (recovery_psi_min <= current_psi <= recovery_psi_max)
      
      if is_near_stable
        # Compute first derivative |ψ'(t)| using central difference
        dt1 = t_history[n] - t_history[n-1]
        psi_deriv = abs((ψ_history[n] - ψ_history[n-1]) / dt1)
        
        # Compute second derivative |ψ''(t)| using central difference
        dt2 = t_history[n-1] - t_history[n-2]
        avg_dt = 0.5 * (dt1 + dt2)
        psi_accel = abs((ψ_history[n] - 2*ψ_history[n-1] + ψ_history[n-2]) / avg_dt^2)
        
        is_recovered = (psi_deriv < recovery_deriv_threshold) && 
                      (psi_accel < recovery_accel_threshold)
        
        if is_recovered
          # println("PERIODIC RECOVERY EVENT @ t=$(round(integrator.t, sigdigits=4)): RECOVERY (A > A*)")
          # println("ψ = $(round(current_psi, sigdigits=4)) ≈ fast.ψ = $(round(states.fast.ψ, sigdigits=4)) (within ±$(recovery_psi_tolerance))")
          # println("|ψ'(t)| = $(round(psi_deriv, sigdigits=3)) < threshold($(recovery_deriv_threshold))")
          # println("|ψ''(t)| = $(round(psi_accel, sigdigits=3)) < threshold($(recovery_accel_threshold))")
          # println("Progress: $(round(100 * (integrator.t - t0)/(t1 - t0), digits=1))% of simulation time")
          
          DiffEqBase.terminate!(integrator)
        end
      end
    end
  end
  
  # Create callbacks
  quench_callback = ContinuousCallback(quench_condition, quench_affect!)
  recovery_period = (t1-t0)/10.0  # Check every 10% of simulation time (reduced frequency for performance)
  recovery_callback = PeriodicCallback(recovery_affect!, recovery_period)
  
  # Combine callbacks
  combined_callback = CallbackSet(quench_callback, recovery_callback)
  
  # Solve with both Sundials events and periodic recovery checks
  sol = solve(prob2, CVODE_BDF(linear_solver = :GMRES);
              abstol = atol, reltol = rtol,
              callback = combined_callback,
              save_everystep = false, save_start = false, save_end = true)
  
  # Compute final ψ
  u_final = @view sol.u[end][1, :]
  ψ_final = Ψ(u_final; ū = states.ū)
  
  # Report final result
  if sol.retcode == :Terminated
    # println("EARLY TERMINATED: Final ψ = $(round(ψ_final, sigdigits=6))")
  else
    # println("RAN TO COMPLETION (near threshold A*): Final ψ = $(round(ψ_final, sigdigits=6))")
    # println(" Target ψ_slow = $(round(states.slow.ψ, sigdigits=6))")
  end
  
  return ψ_final
end

function G(q; f=F, states, lims, fail=+1.0)
  xs, θ   = q
  target  = states.slow.ψ            # the ψ of the UNSTABLE pulse
  prob    = IntervalNonlinearProblem(
    (u,p) -> f((p[1],p[2],u); states=states) - p[3],
    lims,                           
    [xs, θ, target],
  )
  sol = solve(prob, Bisection(); abstol=atol, reltol=rtol)
  if SciMLBase.successful_retcode(sol) && !any(isapprox.(sol.u, lims))
    A = sol.u[1]
  else
    A = fail
  end
  return A
end


end # module
