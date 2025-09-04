#pragma once

#include "fitzhugh_nagumo.h"

class WaveAnalysis {
private:
    Vector xi;
    double xi_min, xi_max;
    
public:
    WaveAnalysis(const Vector& grid);
    
    double compute_psi(const Vector& u, const Vector& u_rest) const;
    double compute_psi(const Matrix& U, const Vector& u_rest) const;
};
