#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <array>
#include <memory>

// Forward declare SUNDIALS types
struct _generic_N_Vector;
typedef struct _generic_N_Vector* N_Vector;

namespace FHN {
    constexpr double Lglob = 2700.0;
    constexpr int N = 1 + (1 << 13);  // 1 + 2^13 = 8193
    constexpr double h = Lglob / (N - 1);
    constexpr double atol = 1e-8;
    constexpr double rtol = 1e-6;
}

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

// Wave structure
struct Wave {
    Matrix u;           // u field (2 x N matrix)
    double psi;
    Vector p;           // parameters [p1, p2, p3]
    
    Wave(const Matrix& u_field, double psi_val, const Vector& params)
        : u(u_field), psi(psi_val), p(params) {}
};

// States structure 
struct States {
    Vector u_rest;     
    Wave fast;          // fast wave
    Wave slow;          // slow wave
    std::pair<double, double> tspan;  // time span for integration
    Matrix u0;          // initial condition matrix
    Vector params;      // PDE parameters
    
    States(const Vector& rest, const Wave& fast_wave, const Wave& slow_wave,
           std::pair<double, double> time_span, const Matrix& initial, const Vector& pde_params)
        : u_rest(rest), fast(fast_wave), slow(slow_wave), tspan(time_span), 
          u0(initial), params(pde_params) {}
};

// Forward declarations
class FiniteDifferenceOperator;
class WaveAnalysis;

// Solver class declaration
class FitzHughNagumoSolver {
private:
    std::unique_ptr<FiniteDifferenceOperator> fd_op;
    std::unique_ptr<WaveAnalysis> wave_analyzer;
    
public:
    FitzHughNagumoSolver();
    ~FitzHughNagumoSolver() = default;
    

    static void fhn_local(double* dx, const double* x, const Vector& p);
    

    static int fhn_pde_rhs(double t, N_Vector y, N_Vector ydot, void* user_data);
    

    Matrix X(const std::vector<double>& q) const;
    

    double F(const std::vector<double>& q, const States& states, bool early_termination = false);
    

    Vector rest(const Vector& p) const;
    
    // Helper functions
    const FiniteDifferenceOperator& get_fd_operator() const { return *fd_op; }
    const WaveAnalysis& get_wave_analyzer() const { return *wave_analyzer; }
};
