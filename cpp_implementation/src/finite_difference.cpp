#include "finite_difference.h"
#include <iostream>

FiniteDifferenceOperator::FiniteDifferenceOperator(int grid_points, double spacing) 
    : N(grid_points), h(spacing) {
    
    // Create spatial grid - range(-Lglob/2, Lglob/2, length=N)
    xi = Vector::LinSpaced(N, -FHN::Lglob/2.0, FHN::Lglob/2.0);
    
    build_laplacian_6th_order_periodic();
}

void FiniteDifferenceOperator::build_laplacian_6th_order_periodic() {
    // CenteredDifference(2, 6, h, N) with PeriodicBC
    // 6th order centered difference stencil for 2nd derivative:
    // Coefficients from finite difference tables
    std::vector<double> coeffs = {
        1.0/90.0,    // i-3
        -3.0/20.0,   // i-2  
        3.0/2.0,     // i-1
        -49.0/18.0,  // i (center)
        3.0/2.0,     // i+1
        -3.0/20.0,   // i+2
        1.0/90.0     // i+3
    };
    
    std::vector<int> offsets = {-3, -2, -1, 0, 1, 2, 3};
    
    std::vector<Eigen::Triplet<double>> triplets;
    double h2_inv = 1.0 / (h * h);
    
    // Build sparse matrix with periodic boundary conditions
    for (int i = 0; i < N; ++i) {
        for (size_t k = 0; k < coeffs.size(); ++k) {
            // Periodic boundary conditions: wrap around using modulo
            int j = (i + offsets[k] + N) % N;
            triplets.push_back(Eigen::Triplet<double>(i, j, coeffs[k] * h2_inv));
        }
    }
    
    laplacian.resize(N, N);
    laplacian.setFromTriplets(triplets.begin(), triplets.end());
}

Vector FiniteDifferenceOperator::apply(const Vector& u) const {
    return laplacian * u;
}
