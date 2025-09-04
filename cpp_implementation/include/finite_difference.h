#pragma once

#include "fitzhugh_nagumo.h"

class FiniteDifferenceOperator {
private:
    int N;
    double h;
    SparseMatrix laplacian;
    Vector xi;  
    
    void build_laplacian_6th_order_periodic();
    
public:
    FiniteDifferenceOperator(int grid_points, double spacing);
    
    const SparseMatrix& get_laplacian() const { return laplacian; }
    const Vector& get_grid() const { return xi; }

    Vector apply(const Vector& u) const;
};
