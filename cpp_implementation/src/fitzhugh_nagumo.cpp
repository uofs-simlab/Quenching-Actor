#include "fitzhugh_nagumo.h"
#include "finite_difference.h"
#include "wave_analysis.h"
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_context.h>
#include <iostream>
#include <cmath>
#include <algorithm>

// User data structure to pass to SUNDIALS
struct UserData {
    SparseMatrix laplacian;
    Vector params;
    int N;
    
    UserData(const SparseMatrix& lap, const Vector& p, int grid_size) 
        : laplacian(lap), params(p), N(grid_size) {}
};

// Early termination data structure
struct EarlyTerminationData {
    double rest_threshold;
    double recovery_psi_tolerance;
    double recovery_deriv_threshold;
    double recovery_accel_threshold;  // Added: second derivative threshold
    double min_time_fraction;
    double min_time;
    double fast_psi;
    Vector u_rest;
    std::vector<double> psi_history;
    std::vector<double> t_history;
    int history_length;
    WaveAnalysis* wave_analyzer;
    
    EarlyTerminationData(double t0, double t1, double fast_psi_val, const Vector& rest_state, WaveAnalysis* analyzer)
        : rest_threshold(0.01), recovery_psi_tolerance(0.5), recovery_deriv_threshold(1e-4),
          recovery_accel_threshold(1e-5), min_time_fraction(0.2), min_time(t0 + min_time_fraction * (t1 - t0)),
          fast_psi(fast_psi_val), u_rest(rest_state), history_length(5), wave_analyzer(analyzer) {}
};

// Data for SUNDIALS with early termination
struct CombinedUserData {
    UserData solver_data;
    EarlyTerminationData* early_term_data;
    
    CombinedUserData(const SparseMatrix& lap, const Vector& p, int grid_size, EarlyTerminationData* et_data)
        : solver_data(lap, p, grid_size), early_term_data(et_data) {}
};

int early_termination_root(double t, N_Vector y, double* gout, void* user_data);

FitzHughNagumoSolver::FitzHughNagumoSolver() {
    fd_op = std::make_unique<FiniteDifferenceOperator>(FHN::N, FHN::h);
    wave_analyzer = std::make_unique<WaveAnalysis>(fd_op->get_grid());
}

void FitzHughNagumoSolver::fhn_local(double* dx, const double* x, const Vector& p) {
    // dx[1] = x[1] * (1 - x[1]) * (x[1] - p[2]) - x[2]
    // dx[2] = p[3] * (p[1] * x[1] - x[2])
    
    dx[0] = x[0] * (1.0 - x[0]) * (x[0] - p(1)) - x[1];  
    dx[1] = p(2) * (p(0) * x[0] - x[1]);
}

int FitzHughNagumoSolver::fhn_pde_rhs(double t, N_Vector y, N_Vector ydot, void* user_data) {
    
    UserData* data;
    CombinedUserData* combined_data = static_cast<CombinedUserData*>(user_data);
    
    // Try to determine if this is CombinedUserData by checking if early_term_data is valid
    // For safety, we'll assume it's always CombinedUserData now
    data = &(combined_data->solver_data);
    
    double* y_data = N_VGetArrayPointer(y);
    double* ydot_data = N_VGetArrayPointer(ydot);
    
    int N = data->N;
    const Vector& p = data->params;
    
    // y is flattened: [u1_1, u1_2, ..., u1_N, u2_1, u2_2, ..., u2_N]
    // First N elements are u component, next N are v component
    
    // Apply reaction terms to each spatial point
    for (int i = 0; i < N; ++i) {
        double x[2] = {y_data[i], y_data[i + N]};  // u[i], v[i]
        double dx[2];
        
        fhn_local(dx, x, p);
        
        ydot_data[i] = dx[0];      // du/dt for point i
        ydot_data[i + N] = dx[1];  // dv/dt for point i
    }
    
    // Add diffusion term: mul!(du, \Delta, x[1, :], 1.0, 1.0)
    // This means: du += \Delta * u (where u is the first component)
    Vector u_vec(N);
    for (int i = 0; i < N; ++i) {
        u_vec(i) = y_data[i];
    }
    
    Vector diffusion = data->laplacian * u_vec;
    
    // Add diffusion to du/dt
    for (int i = 0; i < N; ++i) {
        ydot_data[i] += diffusion(i);
    }
    
    return 0; 
}

Matrix FitzHughNagumoSolver::X(const std::vector<double>& q) const {

    if (q.size() != 3) {
        throw std::runtime_error("X function requires exactly 3 parameters: xs, θ, A");
    }
    
    double xs = q[0], theta = q[1], A = q[2];
    Matrix out = Matrix::Zero(2, FHN::N);
    
    const Vector& xi = fd_op->get_grid();
    
    for (int i = 0; i < FHN::N; ++i) {
        double xi_val = xi(i);
        double arg1 = xi_val - theta + xs / 2.0;
        double arg2 = xi_val - theta - xs / 2.0;
        
        // H = sign function
        double H1 = (arg1 >= 0) ? 1.0 : -1.0;
        double H2 = (arg2 >= 0) ? 1.0 : -1.0;
        
        out(0, i) = (A / 4.0) * (1.0 + H1) * (1.0 - H2);
        // out(1, i) remains zero 
    }
    
    return out;
}

double FitzHughNagumoSolver::F(const std::vector<double>& q, const States& states, bool early_termination) {
    
    // Create perturbation
    Matrix perturbation = X(q);
    
    // Initial condition: u0 + X(q)
    Matrix u0_perturbed = states.u0 + perturbation;
    
    // Create SUNDIALS context (required for SUNDIALS 6.x)
    SUNContext sunctx;
    int flag = SUNContext_Create(NULL, &sunctx);
    if (flag != 0) {
        throw std::runtime_error("SUNContext_Create failed");
    }
    
    // Flatten initial condition for SUNDIALS (u components first, then v components)
    int N = FHN::N;
    N_Vector y0 = N_VNew_Serial(2 * N, sunctx);
    double* y0_data = N_VGetArrayPointer(y0);
    
    for (int i = 0; i < N; ++i) {
        y0_data[i] = u0_perturbed(0, i);      // u component
        y0_data[i + N] = u0_perturbed(1, i);  // v component
    }
    
    // Set up SUNDIALS solver
    void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (!cvode_mem) {
        N_VDestroy(y0);
        SUNContext_Free(&sunctx);
        throw std::runtime_error("CVodeCreate failed");
    }
    
    // Set up data and early termination if enabled
    EarlyTerminationData* et_data = nullptr;
    CombinedUserData* combined_user_data = nullptr;
    
    if (early_termination) {
        // Create early termination data
        et_data = new EarlyTerminationData(states.tspan.first, states.tspan.second, 
                                           states.fast.psi, states.u_rest, wave_analyzer.get());
        
        combined_user_data = new CombinedUserData(fd_op->get_laplacian(), states.params, N, et_data);
    } else {
        combined_user_data = new CombinedUserData(fd_op->get_laplacian(), states.params, N, nullptr);
    }
    
    // Initialize CVODE
    flag = CVodeInit(cvode_mem, fhn_pde_rhs, states.tspan.first, y0);
    if (flag != CV_SUCCESS) {
        if (et_data) delete et_data;
        if (combined_user_data) delete combined_user_data;
        CVodeFree(&cvode_mem);
        N_VDestroy(y0);
        SUNContext_Free(&sunctx);
        throw std::runtime_error("CVodeInit failed");
    }
    
    // Set tolerances
    flag = CVodeSStolerances(cvode_mem, FHN::rtol, FHN::atol);
    if (flag != CV_SUCCESS) {
        if (et_data) delete et_data;
        if (combined_user_data) delete combined_user_data;
        CVodeFree(&cvode_mem);
        N_VDestroy(y0);
        SUNContext_Free(&sunctx);
        throw std::runtime_error("CVodeSStolerances failed");
    }
    
    flag = CVodeSetUserData(cvode_mem, combined_user_data);
    if (flag != CV_SUCCESS) {
        if (et_data) delete et_data;
        if (combined_user_data) delete combined_user_data;
        CVodeFree(&cvode_mem);
        N_VDestroy(y0);
        SUNContext_Free(&sunctx);
        throw std::runtime_error("CVodeSetUserData failed");
    }
    
    // Set up event detection for early termination
    if (early_termination) {
        flag = CVodeRootInit(cvode_mem, 1, early_termination_root);  // 1 root function
        if (flag != CV_SUCCESS) {
            if (et_data) delete et_data;
            if (combined_user_data) delete combined_user_data;
            CVodeFree(&cvode_mem);
            N_VDestroy(y0);
            SUNContext_Free(&sunctx);
            throw std::runtime_error("CVodeRootInit failed");
        }
    }
    
    flag = CVodeSetMaxNumSteps(cvode_mem, 100000); 
    if (flag != CV_SUCCESS) {
        if (et_data) delete et_data;
        if (combined_user_data) delete combined_user_data;
        CVodeFree(&cvode_mem);
        N_VDestroy(y0);
        SUNContext_Free(&sunctx);
        throw std::runtime_error("CVodeSetMaxNumSteps failed");
    }
    
    // Set maximum step size to prevent overly small steps
    flag = CVodeSetMaxStep(cvode_mem, 1.0);  // Max step size of 1.0
    if (flag != CV_SUCCESS) {
        if (et_data) delete et_data;
        if (combined_user_data) delete combined_user_data;
        CVodeFree(&cvode_mem);
        N_VDestroy(y0);
        SUNContext_Free(&sunctx);
        throw std::runtime_error("CVodeSetMaxStep failed");
    }
    
    // Create matrix and linear solver 
    SUNLinearSolver LS = SUNLinSol_SPGMR(y0, PREC_NONE, 0, sunctx);
    
    // Attach linear solver (no matrix needed for SPGMR)
    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (flag != CV_SUCCESS) {
        SUNLinSolFree(LS);
        if (et_data) delete et_data;
        if (combined_user_data) delete combined_user_data;
        CVodeFree(&cvode_mem);
        N_VDestroy(y0);
        SUNContext_Free(&sunctx);
        throw std::runtime_error("CVodeSetLinearSolver failed");
    }
    
    // Solve to final time
    double t_final = states.tspan.second;
    double t_current = states.tspan.first;
    
    flag = CVode(cvode_mem, t_final, y0, &t_current, CV_NORMAL);
    
    // Check if early termination occurred
    bool early_terminated = false;
    if (early_termination && flag == CV_ROOT_RETURN) {
        early_terminated = true;
        // std::cout << "Early termination triggered at t=" << t_current 
        //           << " (" << 100.0 * (t_current - states.tspan.first) / (t_final - states.tspan.first) 
        //           << "% of simulation)" << std::endl;
    } else if (flag < 0) {
        std::cout << "Warning: CVode failed with flag " << flag << std::endl;
    }
    Vector u_final(N);
    for (int i = 0; i < N; ++i) {
        u_final(i) = y0_data[i];  // Extract u component
    }
    
    double final_psi = wave_analyzer->compute_psi(u_final, states.u_rest);
    
    // Cleanup
    SUNLinSolFree(LS);
    CVodeFree(&cvode_mem);
    N_VDestroy(y0);
    SUNContext_Free(&sunctx);
    
    // Clean up early termination data
    if (et_data) delete et_data;
    if (combined_user_data) delete combined_user_data;
    
            // std::cout << "X final value: " << X_final << " for g=" << g << ", theta=" << theta << std::endl;
    
    return final_psi;
}

Vector FitzHughNagumoSolver::rest(const Vector& p) const {
    
    Vector u0(2);
    u0 << 0.0, 0.0;  
    
    // Simple fixed-point iteration 
    for (int iter = 0; iter < 100; ++iter) {
        double dx[2];
        double x[2] = {u0(0), u0(1)};
        
        fhn_local(dx, x, p);
        
        u0(0) -= 0.01 * dx[0];
        u0(1) -= 0.01 * dx[1];
        
        if (std::abs(dx[0]) + std::abs(dx[1]) < 1e-10) {
            break;
        }
    }
    
    return u0;
}

// Root function for SUNDIALS event detection (early termination)
int early_termination_root(double t, N_Vector y, double* gout, void* user_data) {
    CombinedUserData* combined_data = static_cast<CombinedUserData*>(user_data);
    EarlyTerminationData* et_data = combined_data->early_term_data;
    
    if (!et_data) {
        gout[0] = 1.0;  // No event
        return 0;
    }
    
    // Only check after minimum time
    if (t < et_data->min_time) {
        gout[0] = 1.0;  // No event
        return 0;
    }
    
    double* y_data = N_VGetArrayPointer(y);
    int N = combined_data->solver_data.N;
    
    // Extract u component and compute max deviation from rest state
    double max_u_deviation = 0.0;
    for (int i = 0; i < N; ++i) {
        double deviation = std::abs(y_data[i]);  // |u - 0| (rest state is approximately 0)
        max_u_deviation = std::max(max_u_deviation, deviation);
    }
    
    // QUENCHING DETECTION: Simple and fast
    if (max_u_deviation < et_data->rest_threshold) {
        gout[0] = -1.0;  // Trigger event (quenching detected)
        return 0;
    }
    
    // RECOVERY DETECTION: Check if \psi is close to stable wave integral
    Vector u_component(N);
    for (int i = 0; i < N; ++i) {
        u_component(i) = y_data[i];
    }
    
    double current_psi = et_data->wave_analyzer->compute_psi(u_component, et_data->u_rest);
    
    // Store history for derivative computation
    et_data->psi_history.push_back(current_psi);
    et_data->t_history.push_back(t);
    
    // Keep only recent history
    if (et_data->psi_history.size() > static_cast<size_t>(et_data->history_length)) {
        et_data->psi_history.erase(et_data->psi_history.begin());
        et_data->t_history.erase(et_data->t_history.begin());
    }

    // Check if \psi is close to stable wave integral (fast.\psi)
    double recovery_psi_min = et_data->fast_psi - et_data->recovery_psi_tolerance;
    double recovery_psi_max = et_data->fast_psi + et_data->recovery_psi_tolerance;
    bool is_near_stable = (recovery_psi_min <= current_psi && current_psi <= recovery_psi_max);
    
    if (is_near_stable && et_data->psi_history.size() >= 3) {
        size_t n = et_data->psi_history.size();
        
        // Compute first derivative |\psi'(t)| using central difference
        double dt1 = et_data->t_history[n-1] - et_data->t_history[n-2];
        double psi_deriv = std::abs((et_data->psi_history[n-1] - et_data->psi_history[n-2]) / dt1);

        // Compute second derivative |\psi''(t)| using central difference
        double dt2 = et_data->t_history[n-2] - et_data->t_history[n-3];
        double avg_dt = 0.5 * (dt1 + dt2);
        double psi_accel = std::abs((et_data->psi_history[n-1] - 2*et_data->psi_history[n-2] + et_data->psi_history[n-3]) / (avg_dt * avg_dt));

        // RECOVERY DETECTION: Near stable \psi AND small derivatives (wave has stabilized)
        // both first and second derivatives must be small
        bool is_recovered = (psi_deriv < et_data->recovery_deriv_threshold) && 
                           (psi_accel < et_data->recovery_accel_threshold);
        
        if (is_recovered) {
            gout[0] = -1.0;  // Trigger event (recovery detected)
            return 0;
        }
    }
    
    gout[0] = 1.0;  // No event
    return 0;
}
