#include "wave_loader.h"
#include "finite_difference.h"
#include "wave_analysis.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

Matrix WaveLoader::load_wave_data(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::vector<std::vector<double>> data;
    std::string line;
    
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::istringstream iss(line);
        double value;
        
        while (iss >> value) {
            row.push_back(value);
        }
        
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    
    if (data.empty()) {
        throw std::runtime_error("No data found in file: " + filename);
    }
    
    // Convert to Eigen matrix
    int rows = data.size();
    int cols = data[0].size();
    Matrix result(rows, cols);
    
    for (int i = 0; i < rows; ++i) {
        if (data[i].size() != cols) {
            throw std::runtime_error("Inconsistent row sizes in file: " + filename);
        }
        for (int j = 0; j < cols; ++j) {
            result(i, j) = data[i][j];
        }
    }
    
    return result;
}

Vector WaveLoader::load_parameters(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::vector<double> params;
    double value;
    
    while (file >> value) {
        params.push_back(value);
    }
    
    Vector result(params.size());
    for (size_t i = 0; i < params.size(); ++i) {
        result(i) = params[i];
    }
    
    return result;
}

std::tuple<Matrix, Vector, double, Vector> WaveLoader::load_and_center_simple(
    const std::string& Ufile, 
    const std::string& Pfile, 
    const FitzHughNagumoSolver& solver) {
    
    // Load raw data
    Matrix M = load_wave_data(Ufile);
    Vector pf = load_parameters(Pfile);
    
    Vector x = M.col(0);   // spatial coordinate
    Vector u1 = M.col(1);  
    Vector u2 = M.col(2);  
    
    double L = pf(4);  
    
    // Find x0 = x[argmax(u1)] - the position of maximum u1
    int max_idx = 0;
    double max_u1 = u1(0);
    for (int i = 1; i < u1.size(); ++i) {
        if (u1(i) > max_u1) {
            max_u1 = u1(i);
            max_idx = i;
        }
    }
    double x0 = x(max_idx);
    
    // Create U matrix for our grid using proper interpolation
    Matrix U = Matrix::Zero(2, FHN::N);
    
    const auto& fd_op = solver.get_fd_operator();
    const Vector& xi = fd_op.get_grid();

    // Create periodic tiled domain and center at x0
    std::vector<double> x_tiled, u1_tiled, u2_tiled;
    
    // Tile the domain: xtile = vcat(x, x .+ L)
    for (int i = 0; i < x.size(); ++i) {
        x_tiled.push_back(x(i) - x0);  // Center at x0
        u1_tiled.push_back(u1(i));
        u2_tiled.push_back(u2(i));
    }
    for (int i = 0; i < x.size(); ++i) {
        x_tiled.push_back(x(i) + L - x0);  // Periodic extension
        u1_tiled.push_back(u1(i));
        u2_tiled.push_back(u2(i));
    }
    
    // Simple linear interpolation
    for (int i = 0; i < FHN::N; ++i) {
        double xi_val = xi(i);
        
        // Find bracketing points in tiled domain
        int idx = 0;
        while (idx < x_tiled.size() - 1 && x_tiled[idx + 1] < xi_val) {
            idx++;
        }
        
        if (idx >= x_tiled.size() - 1) {
            // Use last value
            U(0, i) = u1_tiled.back();
            U(1, i) = u2_tiled.back();
        } else {
            // Linear interpolation
            double t = (xi_val - x_tiled[idx]) / (x_tiled[idx + 1] - x_tiled[idx]);
            U(0, i) = u1_tiled[idx] + t * (u1_tiled[idx + 1] - u1_tiled[idx]);
            U(1, i) = u2_tiled[idx] + t * (u2_tiled[idx + 1] - u2_tiled[idx]);
        }
    }
    
    Vector p = pf.head(3);
    
    Vector u_rest = solver.rest(p);
    

    const auto& wave_analyzer = solver.get_wave_analyzer();
    Vector u_row = U.row(0);  
    double psi = wave_analyzer.compute_psi(u_row, u_rest);
    
    return std::make_tuple(U, p, psi, u_rest);
}

States WaveLoader::create_states_from_files(
    const std::string& fastU, const std::string& fastP,
    const std::string& slowU, const std::string& slowP,
    const FitzHughNagumoSolver& solver) {
    
    auto [Ufast, pfast, psi_fast, u_rest] = load_and_center_simple(fastU, fastP, solver);
    auto [Uslow, pslow, psi_slow, u_rest_slow] = load_and_center_simple(slowU, slowP, solver);
    
    // Create wave objects
    Wave fast_wave(Ufast, psi_fast, pfast);
    Wave slow_wave(Uslow, psi_slow, pslow);
    
    // Compute time span
    Vector pf_full = load_parameters(fastP);
    double wave_speed = (pf_full.size() > 3) ? pf_full(3) : 1.0;
    std::pair<double, double> tspan = {0.0, FHN::Lglob / (2.0 * wave_speed)};
    
    return States(u_rest, fast_wave, slow_wave, tspan, Ufast, pfast);
}
