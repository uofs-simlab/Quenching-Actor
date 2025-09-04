#include "wave_analysis.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <fstream>
#include <iostream>
#include <stdexcept>

WaveAnalysis::WaveAnalysis(const Vector& grid) : xi(grid) {
    xi_min = xi.minCoeff();
    xi_max = xi.maxCoeff();
}

double WaveAnalysis::compute_psi(const Vector& u, const Vector& u_rest) const {

    
    if (u.size() != xi.size()) {
        throw std::runtime_error("Size mismatch between u and xi in compute_psi");
    }
    
    // Compute integrand: |u - u_rest[0]|
    Vector integrand(u.size());
    for (int i = 0; i < u.size(); ++i) {
        integrand(i) = std::abs(u(i) - u_rest(0));
    }
    
    // Create GSL spline (cubic spline, k=5 in Julia corresponds to high-order spline)
    // GSL doesn't have exact 5th order, so we use cubic which is robust
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, u.size());
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    
    if (!spline || !acc) {
        throw std::runtime_error("Failed to allocate GSL spline or accelerator");
    }
    
    // Initialize spline with our data
    int status = gsl_spline_init(spline, xi.data(), integrand.data(), u.size());
    if (status != GSL_SUCCESS) {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
        throw std::runtime_error("Failed to initialize GSL spline");
    }
    
    // Integrate spline over entire domain
    double result = gsl_spline_eval_integ(spline, xi_min, xi_max, acc);
    
    // Cleanup
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    
    return result;
}

double WaveAnalysis::compute_psi(const Matrix& U, const Vector& u_rest) const {

    double total_psi = 0.0;
    
    for (int row = 0; row < U.rows(); ++row) {
        Vector u_row = U.row(row);
        total_psi += compute_psi(u_row, u_rest);
    }
    
    return total_psi;
}
