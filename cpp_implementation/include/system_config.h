#ifndef SYSTEM_CONFIG_H
#define SYSTEM_CONFIG_H

#include "job_structures.h"
#include "bracket_optimizer.h"

struct system_config {
    // Computation options (set from command-line config, not hardcoded here)
    bool enable_bracket_optimization = false;  // Will be set from config.enable_bracket
    bool enable_early_termination = false;     // Will be set from config.enable_early_termination
    
    // Grid parameters
    float theta_step = 10.0f;       // Dynamic: 10.0f, Static: 20.0f
    float xs_step = 50.0f;          // Both: 50.0f (calculated)
    float max_theta = 40.0f;        // Dynamic: 40.0f, Static: 80.0f
    float xs_min = 0.1f;            // Dynamic: 0.1f, Static: 0.05f
    float xs_max = 700.0f;          // Dynamic: 700.0f, Static: 1000.0f
    int NDNS = 10;                  // Both: 10
    
    // Refinement parameters
    int max_refinement_depth = 2;   // Dynamic: 2, Static: 8
    
    // Actor system parameters
    int actor_number = 31;          // Both: 31
    
    // Bracket optimization parameters
    bracket_optimizer::BracketConfig bracket_config;
    
    // System type
    enum SystemType { DYNAMIC, STATIC } system_type = DYNAMIC;
    
    // Factory methods for creating standard configurations
    static system_config create_dynamic_config(bool enable_bracket = false, bool enable_early_term = false) {
        system_config config;
        config.system_type = DYNAMIC;
        config.enable_bracket_optimization = enable_bracket;
        config.enable_early_termination = enable_early_term;
        config.theta_step = 10.0f;
        config.max_theta = 40.0f;
        config.xs_min = 0.1f;
        config.xs_max = 700.0f;
        config.max_refinement_depth = 2;
        return config;
    }
    
    static system_config create_static_config(bool enable_bracket = false, bool enable_early_term = false) {
        system_config config;
        config.system_type = STATIC;
        config.enable_bracket_optimization = enable_bracket;
        config.enable_early_termination = enable_early_term;
        config.theta_step = 20.0f;
        config.max_theta = 80.0f;
        config.xs_min = 0.05f;
        config.xs_max = 1000.0f;
        config.max_refinement_depth = 8;
        return config;
    }
    
    // Calculate derived parameters
    void finalize() {
        xs_step = (xs_max - xs_min) / (NDNS - 1);
    }
};

#endif // SYSTEM_CONFIG_H
