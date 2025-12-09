#include "fitzhugh_nagumo.h"
#include "finite_difference.h"
#include "wave_analysis.h"
#include "wave_loader.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_context.h>

// Function to load wave data from file (mirrors Julia's readdlm)
Matrix load_wave_data(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::vector<std::vector<double>> data;
    std::string line;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> row;
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
    
    int rows = data.size();
    int cols = data[0].size();
    Matrix result(rows, cols);
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result(i, j) = data[i][j];
        }
    }
    
    return result;
}

// Function to load parameter vector from file
Vector load_parameters(const std::string& filename) {
    Matrix param_matrix = load_wave_data(filename);
    
    Vector params;
    if (param_matrix.rows() > param_matrix.cols()) {
        params = param_matrix.col(0);
    } else {
        params = param_matrix.row(0);
    }
    
    return params;
}

std::tuple<Matrix, Vector, double, Vector> load_and_center_simple(const std::string& Ufile, const std::string& Pfile, const FitzHughNagumoSolver& solver) {
    std::cout << "Loading files: " << Ufile << ", " << Pfile << std::endl;
    
    // Use proper WaveLoader for correct centering
    WaveLoader loader;
    std::cout << "WaveLoader created, calling load_and_center_simple..." << std::endl;
    auto result = loader.load_and_center_simple(Ufile, Pfile, solver);
    std::cout << "WaveLoader completed successfully" << std::endl;
    return result;
}

States create_states_from_files(const std::string& fastU, const std::string& fastP,
                               const std::string& slowU, const std::string& slowP,
                               const FitzHughNagumoSolver& solver) {
    auto [Ufast, pfast, psi_fast, u_rest] = load_and_center_simple(fastU, fastP, solver);
    auto [Uslow, pslow, psi_slow, u_rest2] = load_and_center_simple(slowU, slowP, solver);
    
    std::cout << "Parameter sizes: pfast=" << pfast.size() << ", pslow=" << pslow.size() << std::endl;
    if (pfast.size() > 0) std::cout << "pfast[0-2]: " << pfast(0) << ", " << pfast(1) << ", " << pfast(2) << std::endl;
    if (pfast.size() > 3) std::cout << "pfast[3]: " << pfast(3) << std::endl;
    
    Wave fast_wave(Ufast, psi_fast, pfast);
    Wave slow_wave(Uslow, psi_slow, pslow);
    
    // Time span calculation: tspan = (0.0, Lglob / (2 * pfast[4]))
    double Lglob = FHN::Lglob;
    double p4_value = (pfast.size() > 3) ? pfast(3) : 1.0; // Default if not available
    std::pair<double, double> tspan = {0.0, Lglob / (2.0 * p4_value)}; // pfast[4] in Julia is pfast(3) in C++ (0-indexed)
    std::cout << "Using p4_value=" << p4_value << " for tspan calculation" << std::endl;
    
    // Debug: print dimensions
    std::cout << "Creating States with Ufast dimensions: " << Ufast.rows() << " x " << Ufast.cols() << std::endl;
    std::cout << "Expected grid size N = " << FHN::N << std::endl;
    
    return States(u_rest, fast_wave, slow_wave, tspan, Ufast, pfast);
}

// Modified PDE solver that saves snapshots
class SnapshotSolver {
private:
    FitzHughNagumoSolver base_solver;
    std::vector<Matrix> snapshots;
    std::vector<double> time_points;
    int snapshot_count;
    int max_snapshots;
    
public:
    SnapshotSolver(int max_snaps = 100) : max_snapshots(max_snaps), snapshot_count(0) {
        snapshots.reserve(max_snapshots);
        time_points.reserve(max_snapshots);
    }
    
    // Callback function for taking snapshots
    static int snapshot_callback(double t, N_Vector y, N_Vector ydot, void* user_data) {
        SnapshotSolver* solver = static_cast<SnapshotSolver*>(user_data);
        return solver->take_snapshot(t, y);
    }
    
    int take_snapshot(double t, N_Vector y) {
        if (snapshot_count >= max_snapshots) {
            return 0; // Don't take more snapshots
        }
        
        double* y_data = N_VGetArrayPointer(y);
        int N = FHN::N;
        
        Matrix snapshot(2, N);
        for (int i = 0; i < N; ++i) {
            snapshot(0, i) = y_data[i];      // u component
            snapshot(1, i) = y_data[i + N];  // v component
        }
        
        snapshots.push_back(snapshot);
        time_points.push_back(t);
        snapshot_count++;
        
        std::cout << "Snapshot " << snapshot_count << " taken at t=" << t << std::endl;
        return 0;
    }
    
    void solve_with_snapshots(const std::vector<double>& q, const States& states) {
        snapshots.clear();
        time_points.clear();
        snapshot_count = 0;
        
        // Create perturbation
        Matrix perturbation = base_solver.X(q);
        
        // Debug: print dimensions and perturbation info
        std::cout << "states.u0 dimensions: " << states.u0.rows() << " x " << states.u0.cols() << std::endl;
        std::cout << "perturbation dimensions: " << perturbation.rows() << " x " << perturbation.cols() << std::endl;
        
        // Debug: Check where perturbation is non-zero
        Vector pert_u = perturbation.row(0);
        double pert_max = pert_u.maxCoeff();
        double pert_min = pert_u.minCoeff();
        int non_zero_count = 0;
        int first_nonzero = -1, last_nonzero = -1;
        for (int i = 0; i < pert_u.size(); ++i) {
            if (std::abs(pert_u(i)) > 1e-12) {
                if (first_nonzero == -1) first_nonzero = i;
                last_nonzero = i;
                non_zero_count++;
            }
        }
        
        const auto& fd_op = base_solver.get_fd_operator();
        const Vector& xi = fd_op.get_grid();
        
        std::cout << "DEBUG PERTURBATION: xs=" << q[0] << ", theta=" << q[1] << ", A=" << q[2] << std::endl;
        std::cout << "DEBUG PERTURBATION: max=" << pert_max << ", min=" << pert_min << std::endl;
        std::cout << "DEBUG PERTURBATION: non-zero points=" << non_zero_count << std::endl;
        if (first_nonzero >= 0) {
            std::cout << "DEBUG PERTURBATION: range xi=[" << xi(first_nonzero) << ", " << xi(last_nonzero) << "]" << std::endl;
            std::cout << "DEBUG PERTURBATION: center should be at xi=" << q[1] << std::endl;
            std::cout << "DEBUG PERTURBATION: width should be xs=" << q[0] << std::endl;
        }
        
        // Ensure dimensions match
        if (states.u0.rows() != perturbation.rows() || states.u0.cols() != perturbation.cols()) {
            throw std::runtime_error("Dimension mismatch between u0 and perturbation");
        }
        
        Matrix u0_perturbed = states.u0 + perturbation;
        
        // Create SUNDIALS context
        SUNContext sunctx;
        int flag = SUNContext_Create(NULL, &sunctx);
        if (flag != 0) {
            throw std::runtime_error("SUNContext_Create failed");
        }
        
        int N = FHN::N;
        N_Vector y0 = N_VNew_Serial(2 * N, sunctx);
        double* y0_data = N_VGetArrayPointer(y0);
        
        // Set initial conditions
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
        
        // User data for PDE right-hand side
        struct UserData {
            SparseMatrix laplacian;
            Vector params;
            int N;
            UserData(const SparseMatrix& lap, const Vector& p, int grid_size)
                : laplacian(lap), params(p), N(grid_size) {}
        };
        
        UserData user_data(base_solver.get_fd_operator().get_laplacian(), states.params, N);
        
        // Initialize CVODE
        flag = CVodeInit(cvode_mem, FitzHughNagumoSolver::fhn_pde_rhs, states.tspan.first, y0);
        if (flag != CV_SUCCESS) {
            CVodeFree(&cvode_mem);
            N_VDestroy(y0);
            SUNContext_Free(&sunctx);
            throw std::runtime_error("CVodeInit failed");
        }
        
        flag = CVodeSStolerances(cvode_mem, FHN::rtol, FHN::atol);
        flag = CVodeSetUserData(cvode_mem, &user_data);
        
        // Set up linear solver
        SUNLinearSolver LS = SUNLinSol_SPGMR(y0, SUN_PREC_NONE, 0, sunctx);
        flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
        
        // Solve with snapshots
        double t_final = states.tspan.second;
        double t_current = states.tspan.first;
        double dt = (t_final - states.tspan.first) / max_snapshots;
        
        // Take initial snapshot
        take_snapshot(t_current, y0);
        
        // Solve step by step to capture snapshots
        for (int i = 1; i < max_snapshots; ++i) {
            double t_next = states.tspan.first + i * dt;
            flag = CVode(cvode_mem, t_next, y0, &t_current, CV_NORMAL);
            
            if (flag < 0) {
                std::cout << "CVode failed at t=" << t_current << " with flag=" << flag << std::endl;
                break;
            }
            
            take_snapshot(t_current, y0);
        }
        
        // Cleanup
        SUNLinSolFree(LS);
        CVodeFree(&cvode_mem);
        N_VDestroy(y0);
        SUNContext_Free(&sunctx);
    }
    
    void save_snapshots(const std::string& filename_prefix) {
        for (int i = 0; i < snapshots.size(); ++i) {
            std::string filename = filename_prefix + "_frame_" + std::to_string(i) + ".dat";
            std::ofstream file(filename);
            
            file << std::scientific << std::setprecision(15);
            file << "# Time: " << time_points[i] << std::endl;
            file << "# Columns: xi, u, v" << std::endl;
            
            const Vector& xi = base_solver.get_fd_operator().get_grid();
            for (int j = 0; j < FHN::N; ++j) {
                file << xi(j) << " " << snapshots[i](0, j) << " " << snapshots[i](1, j) << std::endl;
            }
            
            file.close();
        }
        
        std::cout << "Saved " << snapshots.size() << " snapshots with prefix: " << filename_prefix << std::endl;
    }
    
    const std::vector<Matrix>& get_snapshots() const { return snapshots; }
    const std::vector<double>& get_time_points() const { return time_points; }
};

int main(int argc, char* argv[]) {
    try {
        std::cout << std::fixed << std::setprecision(15);
        std::cout << "=== FHN SNAPSHOT GENERATOR ===" << std::endl;
        
        // Load states
        std::string base = "/globalhome/tus210/HPC/quenchin_actor/";
        FitzHughNagumoSolver solver;
        States states = create_states_from_files(
            base + "waves/index_11/2/U", base + "waves/index_11/2/p",
            base + "waves/index_10/2/U", base + "waves/index_10/2/p",
            solver
        );
        
        std::cout << "Loaded states successfully" << std::endl;
        std::cout << "Time span: [" << states.tspan.first << ", " << states.tspan.second << "]" << std::endl;
        
        // Parameters as specified: gg=2, theta=20, xs=30
        double gg = 2.0;
        double theta = 20.0;
        double xs = 30.0;
        
        // Test with three different perturbations: 0, -0.58, 0.7
        std::vector<double> perturbations = {0.0, -0.58, 0.7};
        
        for (double A : perturbations) {
            std::cout << "\n=== Running simulation with perturbation A = " << A << " ===" << std::endl;
            
            std::vector<double> q = {xs, theta, A};
            std::cout << "Parameters: xs=" << xs << ", θ=" << theta << ", A=" << A << std::endl;
            
            SnapshotSolver snapshot_solver(100);
            
            auto start_time = std::chrono::high_resolution_clock::now();
            snapshot_solver.solve_with_snapshots(q, states);
            auto end_time = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            std::cout << "Simulation completed in " << duration.count() << " ms" << std::endl;
            
            // Save snapshots
            std::string prefix = "fhn_snapshots_A" + std::to_string(A).substr(0, 6);
            std::replace(prefix.begin(), prefix.end(), '.', 'p');
            std::replace(prefix.begin(), prefix.end(), '-', 'n');
            
            snapshot_solver.save_snapshots(prefix);
        }
        
        std::cout << "\n=== ALL SIMULATIONS COMPLETED ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}