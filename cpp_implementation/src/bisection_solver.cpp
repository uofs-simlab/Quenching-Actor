#include "bisection_solver.h"
#include "wave_loader.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>

BisectionSolver::BisectionSolver(FitzHughNagumoSolver& solver, const BisectionParams& params)
    : solver_(solver), params_(params), last_iterations_(0), last_convergence_error_(0.0),
      last_successful_bracket_({0.0, 0.0}) {
}

double BisectionSolver::evaluate_objective(const std::vector<double>& q, const States& states) {
    if (q.size() != 3) {
        throw std::invalid_argument("evaluate_objective: q must have 3 elements [xs, theta, A]");
    }
    
    // Call F function with early termination flag and subtract target (slow.psi)
    double f_result = solver_.F(q, states, params_.enable_early_termination);
    double target = states.slow.psi;
    
    return f_result - target;
}

std::pair<bool, std::pair<double, double>> BisectionSolver::check_sign_change(
    const std::vector<double>& q,
    const States& states,
    double a_min,
    double a_max) {
    
    if (q.size() != 2) {
        throw std::invalid_argument("check_sign_change: q must have 2 elements [xs, theta]");
    }
    
    // Evaluate at boundaries
    std::vector<double> q_min = {q[0], q[1], a_min};
    std::vector<double> q_max = {q[0], q[1], a_max};
    
    double f_min = evaluate_objective(q_min, states);
    double f_max = evaluate_objective(q_max, states);
    
    bool has_sign_change = (f_min * f_max) < 0;
    
    return std::make_pair(has_sign_change, std::make_pair(f_min, f_max));
}

double BisectionSolver::solve(const std::vector<double>& q, 
                             const States& states,
                             double a_min, 
                             double a_max,
                             const std::map<double, double>* precomputed) {
    
    if (q.size() != 2) {
        throw std::invalid_argument("solve: q must have 2 elements [xs, theta]");
    }
    
    double xs = q[0];
    double theta = q[1];
    double target = states.slow.psi;
    
    // Initial function evaluations - use precomputed values if available
    double f_min, f_max;
    
    if (precomputed && precomputed->find(a_min) != precomputed->end()) {
        f_min = precomputed->at(a_min) - target;
    } else {
        std::vector<double> q_min = {xs, theta, a_min};
        f_min = evaluate_objective(q_min, states);
    }
    
    if (precomputed && precomputed->find(a_max) != precomputed->end()) {
        f_max = precomputed->at(a_max) - target;
    } else {
        std::vector<double> q_max = {xs, theta, a_max};
        f_max = evaluate_objective(q_max, states);
    }
    
    // Check for sign change
    if (f_min * f_max > 0) {
        if (params_.verbose) {
            std::cout << "No sign change: f_min=" << f_min << ", f_max=" << f_max << std::endl;
        }
        return params_.fail_value;
    }
    
    // Ensure proper ordering (f_min < 0, f_max > 0)
    if (f_min > 0) {
        std::swap(a_min, a_max);
        std::swap(f_min, f_max);
        if (params_.verbose) {
            std::cout << "Swapped bracket ordering: new bracket=[" << a_min << ", " << a_max << "]" << std::endl;
        }
    }
    
    // Store original limits for boundary check
    double orig_a_min = std::min(a_min, a_max);
    double orig_a_max = std::max(a_min, a_max);
    
    // Bisection loop
    last_iterations_ = 0;
    for (int iteration = 1; iteration <= params_.max_iter; ++iteration) {
        last_iterations_ = iteration;
        
        double a_mid = (a_min + a_max) / 2.0;
        std::vector<double> q_mid = {xs, theta, a_mid};
        double f_mid = evaluate_objective(q_mid, states);
        
        if (params_.verbose) {
            std::cout << "Iteration " << iteration << ": a_mid=" << a_mid 
                      << ", f_mid=" << f_mid << ", bracket_width=" << (a_max - a_min) << std::endl;
        }
        
        // Check convergence
        if (std::abs(f_mid) < params_.tolerance || std::abs(a_max - a_min) < params_.tolerance) {
            // Check if solution is at boundary (treat as failure)
            if (std::abs(a_mid - orig_a_min) < 1e-12 || std::abs(a_mid - orig_a_max) < 1e-12) {
                if (params_.verbose) {
                    std::cout << "Solution at boundary, treating as failure" << std::endl;
                }
                return params_.fail_value;
            }
            
            last_convergence_error_ = std::abs(f_mid);
            last_successful_bracket_ = {a_min, a_max};
            
            if (params_.verbose) {
                std::cout << "Converged: A=" << a_mid << ", f_mid=" << f_mid 
                          << ", iterations=" << iteration << std::endl;
            }
            
            return a_mid;
        }
        
        // Update bracket
        if (f_mid < 0) {
            a_min = a_mid;
            f_min = f_mid;
        } else {
            a_max = a_mid;
            f_max = f_mid;
        }
    }
    
    // Max iterations reached - check final result
    double a_final = (a_min + a_max) / 2.0;
    if (std::abs(a_final - orig_a_min) < 1e-12 || std::abs(a_final - orig_a_max) < 1e-12) {
        if (params_.verbose) {
            std::cout << "Max iterations reached, solution at boundary" << std::endl;
        }
        return params_.fail_value;
    }
    
    last_convergence_error_ = std::abs(a_max - a_min);
    last_successful_bracket_ = {a_min, a_max};
    
    if (params_.verbose) {
        std::cout << "Max iterations reached: A=" << a_final 
                  << ", bracket_width=" << (a_max - a_min) << std::endl;
    }
    
    return a_final;
}

double BisectionSolver::sweep(const std::vector<double>& q,
                             const States& states,
                             const std::vector<std::pair<double, double>>& brackets,
                             const std::map<double, double>* precomputed) {
    
    if (params_.verbose) {
        std::cout << "Starting bisection sweep with " << brackets.size() << " bracket(s)..." << std::endl;
    }
    
    for (size_t i = 0; i < brackets.size(); ++i) {
        double a_min = brackets[i].first;
        double a_max = brackets[i].second;
        
        if (params_.verbose) {
            std::cout << "Trying bracket " << (i+1) << ": [" << a_min << ", " << a_max << "]" << std::endl;
        }
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        double A_try = solve(q, states, a_min, a_max, precomputed);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        if (A_try != params_.fail_value) {
            if (params_.verbose) {
                std::cout << "SUCCESS in bracket " << (i+1) << " after " 
                          << duration.count() << "ms: A* = " << A_try 
                          << " (iterations: " << last_iterations_ << ")" << std::endl;
            }
            return A_try;
        } else {
            if (params_.verbose) {
                std::cout << "FAILED bracket " << (i+1) << " after " 
                          << duration.count() << "ms" << std::endl;
            }
        }
    }
    
    if (params_.verbose) {
        std::cout << "BISECTION FAILED - no valid solution found in any bracket" << std::endl;
    }
    
    return params_.fail_value;
}

double sweep_main(const std::string& fast_U_file,
                    const std::string& fast_P_file,
                    const std::string& slow_U_file,
                    const std::string& slow_P_file,
                    int n,
                    double usmax,
                    double gg,
                    double xs,
                    double theta,
                    bool verbose,
                    bool enable_early_termination) {
    
    // Create solver
    FitzHughNagumoSolver solver;
    
    // Set up bisection parameters
    BisectionSolver::BisectionParams params;
    params.verbose = verbose;
    params.max_iter = 50;
    params.tolerance = 1e-5;
    params.fail_value = 1.0;
    params.enable_early_termination = enable_early_termination;  // Pass early termination flag
    
    BisectionSolver bisection_solver(solver, params);
    
    if (verbose) {
        std::cout << "Parameters: xs=" << xs << ", θ=" << theta << ", n=" << n << std::endl;
    }
    
    // Compute usmin from n: usmin = -1000/(2^n)
    double usmin = -1000.0 / std::pow(2.0, n);
    
    if (verbose) {
        std::cout << "Computed usmin = " << usmin << " (from n=" << n << ")" << std::endl;
    }
    
    // Load data into States object using WaveLoader
    States states = WaveLoader::create_states_from_files(fast_U_file, fast_P_file, slow_U_file, slow_P_file, solver);
    
    std::vector<std::pair<double, double>> brackets;
    
    if (n == 0) {
        // No bracket optimization - use single wide bracket
        brackets.push_back({-1000.0, 0.0});
    } else {
        // Bracket optimization enabled - use multiple brackets
        brackets.push_back({usmin, usmax});      // Optimized bracket
        brackets.push_back({-1000.0, usmin});    // Fallback bracket
    }
    
    // Pre-check for sign change in widest bracket
    std::vector<double> q = {xs, theta};
    
    auto start_time = std::chrono::high_resolution_clock::now();
    auto [has_sign_change, f_values] = bisection_solver.check_sign_change(q, states, -1000.0, 0.0);
    auto end_time = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    if (!has_sign_change) {
        if (verbose) {
            std::cout << "No sign change in bracket - returning failure" << std::endl;
            std::cout << "f(-1000) = " << f_values.first << ", f(0) = " << f_values.second << std::endl;
        }
        return 1.0;
    }
    
    if (verbose) {
        std::cout << "Sign change detected in " << duration.count() << "ms" << std::endl;
        std::cout << "f(-1000) = " << f_values.first << ", f(0) = " << f_values.second 
                  << ", target = " << states.slow.psi << std::endl;
    }
    
    // Store precomputed values to avoid recomputation
    std::map<double, double> precomputed_vals;
    precomputed_vals[-1000.0] = f_values.first + states.slow.psi;  // Add target back
    precomputed_vals[0.0] = f_values.second + states.slow.psi;
    
    double A = bisection_solver.sweep(q, states, brackets, &precomputed_vals);
    
    if (verbose) {
        if (A != 1.0) {
            std::cout << "FINAL RESULT: A* = " << A << std::endl;
            std::cout << "Convergence: " << bisection_solver.get_last_iterations() 
                      << " iterations, error = " << bisection_solver.get_last_convergence_error() << std::endl;
        } else {
            std::cout << "BISECTION FAILED - no valid solution found" << std::endl;
        }
    }
    
    return A;
}
