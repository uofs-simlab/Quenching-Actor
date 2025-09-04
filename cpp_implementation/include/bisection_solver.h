#ifndef BISECTION_SOLVER_H
#define BISECTION_SOLVER_H

#include "fitzhugh_nagumo.h"
#include "finite_difference.h"
#include "wave_analysis.h"
#include <vector>
#include <map>
#include <utility>

/**
 * @brief Bisection solver for finding quenching threshold A
 * 
 * This class implements the bisection algorithm to find the value of A where
 * the FitzHugh-Nagumo system reaches a target \psi value.
 */
class BisectionSolver {
public:
    struct BisectionParams {
        int max_iter;
        double tolerance;
        double fail_value;
        bool verbose;
        bool enable_early_termination;
        
        BisectionParams() : max_iter(50), tolerance(1e-5), fail_value(1.0), verbose(false), enable_early_termination(false) {}
    };
    
    BisectionSolver(FitzHughNagumoSolver& solver, const BisectionParams& params = BisectionParams());
    
    double solve(const std::vector<double>& q, 
                 const States& states,
                 double a_min, 
                 double a_max,
                 const std::map<double, double>* precomputed = nullptr);
                 
    double sweep(const std::vector<double>& q,
                 const States& states,
                 const std::vector<std::pair<double, double>>& brackets,
                 const std::map<double, double>* precomputed = nullptr);
    
    std::pair<bool, std::pair<double, double>> check_sign_change(
        const std::vector<double>& q,
        const States& states,
        double a_min,
        double a_max);
    
    // Getters for diagnostics
    int get_last_iterations() const { return last_iterations_; }
    double get_last_convergence_error() const { return last_convergence_error_; }
    std::pair<double, double> get_last_bracket() const { return last_successful_bracket_; }
    
private:
    FitzHughNagumoSolver& solver_;
    BisectionParams params_;
    
    // Diagnostic info
    int last_iterations_;
    double last_convergence_error_;
    std::pair<double, double> last_successful_bracket_;
    
    double evaluate_objective(const std::vector<double>& q, const States& states);
};


double sweep_main(const std::string& fast_U_file,
                    const std::string& fast_P_file,
                    const std::string& slow_U_file,
                    const std::string& slow_P_file,
                    int n,
                    double usmax,
                    double gg,
                    double xs,
                    double theta,
                    bool verbose = false,
                    bool enable_early_termination = false);

#endif // BISECTION_SOLVER_H
