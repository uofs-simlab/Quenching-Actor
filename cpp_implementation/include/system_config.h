#ifndef SYSTEM_CONFIG_H
#define SYSTEM_CONFIG_H

#include "job_structures.h"
#include "bracket_optimizer.h"

struct system_config {
    // Grid parameters
    float L = 2700.0f;              // System length parameter (same as FHN::Lglob)
    float h = L / (8192.0f);        // Grid spacing: h = Lglob/(N-1) where N = 1 + 2^13 = 8193
    float theta_step = L/4.0f;      // Dynamic: L/4.0f ≈ 675, Static: 20.0f
    float xs_step = 50.0f;          // Placeholder - calculated in finalize() as (xs_max-xs_min)/(NDNS-1)
    float max_theta = L/2.0f;       // Dynamic: L/2.0f ≈ 1350, Static: 80.0f  
    float xs_min = 2*h;             // Dynamic: 2*h ≈ 0.66, Static: 0.05f
    float xs_max = L;               // Dynamic: L ≈ 2700, Static: 1000.0f
    int NDNS = 5;                   // Dynamic: 5, Static: 10
    
    // Refinement parameters
    int max_refinement_level = 8;   // Maximum global refinement levels
    float min_refinement_radius_factor = 50.0f;  
    int connectivity_range = 4;     // Neighbor connectivity range
    float radius_reduction_factor = 0.6f;  // Radius reduction per level (60% retention)
    
    // Base radius calculation parameters
    float base_radius_multiplier = 2.5f;   // Multiplier for minimum spacing to get base radius
    float coarse_grid_radius_multiplier = 1.0f;  // Additional multiplier for coarse grid connectivity
    
    // Actor system parameters
    int actor_number = 31;          // Number of worker actors
    
    // Spatial hash parameters
    float theta_cell_size_factor = 5.0f;   // theta_step / factor
    float xs_cell_size_factor = 3.5f;      // xs_step / factor
    
    // Problem-specific parameters
    std::vector<int> gg_levels = {3, 2, 1, 0};  // gg levels to process  
    
    // Bracket optimization parameters
    bracket_optimizer::BracketConfig bracket_config;
    
    // System type
    enum SystemType { DYNAMIC, STATIC } system_type = DYNAMIC;
    
    // Factory methods for creating standard configurations
    static system_config create_dynamic_config(bool enable_bracket = false, bool enable_early_term = false) {
        system_config conf;
        conf.system_type = DYNAMIC;
        conf.L = 2700.0f;
        conf.h = conf.L / 8192.0f;
        conf.theta_step = conf.L / 4.0f;  // ~675 degrees per step (5 points total)
        conf.max_theta = conf.L / 2.0f;   // Will be 1350.0 (range: [-1350, 1350])
        conf.xs_min = 2 * conf.h;        // Start from 2h ≈ 0.66
        conf.xs_max = conf.L;             // Will be 2700.0
        conf.NDNS = 5;                      // Reduced for much coarser initial grid
        conf.max_refinement_level = 8;     // Maximum global refinement levels
        conf.min_refinement_radius_factor = 50.0f;  // Same as static for consistency
        conf.connectivity_range = 4;       // 4 grid spacings connectivity
        conf.radius_reduction_factor = 0.6f; // 60% retention per level (same as static)
        conf.base_radius_multiplier = 2.5f; // Base radius calculation
        conf.coarse_grid_radius_multiplier = 1.0f; // Coarse grid gets base radius
        conf.actor_number = 31;            // Number of workers
        conf.gg_levels = {3, 2, 1, 0};     // gg levels to process
        return conf;
    }
    
    static system_config create_static_config(bool enable_bracket = false, bool enable_early_term = false) {
        system_config conf;
        conf.system_type = STATIC;
        
        // Static-specific grid parameters
        conf.L = 2700.0f;                  
        conf.h = conf.L / 8192.0f;          
        conf.theta_step = conf.L / 4.0f;            
        conf.max_theta = conf.L / 2.0f;         
        conf.xs_min = 2 * conf.h;               
        conf.xs_max = conf.L;             
        conf.NDNS = 5;                    
        
        // Static-specific refinement parameters
        conf.max_refinement_level = 8;
        conf.min_refinement_radius_factor = 50.0f;  // Same as before
        conf.connectivity_range = 3;
        conf.radius_reduction_factor = 0.6f; // 60% retention per level
        conf.base_radius_multiplier = 2.5f; // Base radius calculation
        conf.coarse_grid_radius_multiplier = 1.0f; // Coarse grid gets base radius
        conf.actor_number = 31;
        conf.gg_levels = {3, 2, 1, 0};
        return conf;
    }
    
    // Calculate derived parameters
    void finalize() {
        xs_step = (xs_max - xs_min) / (NDNS - 1);
    }
};

#endif // SYSTEM_CONFIG_H
