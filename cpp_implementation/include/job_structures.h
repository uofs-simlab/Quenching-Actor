#ifndef JOB_STRUCTURES_H
#define JOB_STRUCTURES_H

#include <string>

/**
 * @brief Job structure containing all parameters needed for a quenching computation
 */
struct job_t {
    int gg;           // Gamma parameter (0-3)
    float theta;      // Angle parameter (-40 to +40 degrees)
    float xs;         // Spatial parameter (0.1 to 700)
    int n;            // usmin/usmax control (0 = -1000, 1 = 0)
    std::string Ufile;  // Fast wave U file path
    std::string Pfile;  // Fast wave P file path
    std::string ufile;  // Slow wave U file path
    std::string pfile;  // Slow wave P file path

    job_t() = default;
    
    job_t(int g, float t, float x, int n_val, 
          const std::string& U, const std::string& P,
          const std::string& u, const std::string& p)
        : gg(g), theta(t), xs(x), n(n_val), Ufile(U), Pfile(P), ufile(u), pfile(p) {}
};

/**
 * @brief Result entry for caching uq values with corresponding xs values
 */
struct result_entry {
    float xs;
    double uq;
    
    result_entry() = default;
    result_entry(float x, double u) : xs(x), uq(u) {}
};

/**
 * @brief Configuration structure for computation options
 */
struct computation_config {
    bool enable_bracket_optimization = true;   // Enable bracket initialization optimization
    bool enable_early_termination = true;      // Enable early termination in solver
    double bracket_trust_margin = 1.5;         // Trust margin for bracket optimization
    int max_refinement_depth = 2;              // Maximum refinement depth for adaptive mesh
    
    computation_config() = default;
};

#endif // JOB_STRUCTURES_H
