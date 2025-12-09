#include "config.h"
#include "system_config.h"
#include "tuple_hash.h" // Hash specializations for tuples
#include "job_structures.h"
#include "bracket_optimizer.h"
#include "fitzhugh_nagumo.h"
#include "bisection_solver.h"
#include "wave_loader.h"
#include <chrono>
#include <deque>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

//JOB & CACHE TYPES - Using common structures

// CAF serialization support for job_t
template <class Inspector>
bool inspect(Inspector& f, job_t& x) {
    return f.object(x).fields(
        f.field("gg", x.gg),
        f.field("theta", x.theta),
        f.field("xs", x.xs),
        f.field("n", x.n),
        f.field("Ufile", x.Ufile),
        f.field("Pfile", x.Pfile),
        f.field("ufile", x.ufile),
        f.field("pfile", x.pfile)
    );
}

CAF_BEGIN_TYPE_ID_BLOCK(my_project, caf::first_custom_type_id)
    CAF_ADD_TYPE_ID(my_project, (job_t))
CAF_END_TYPE_ID_BLOCK(my_project)

/*
==================================================================================================================
SPATIAL HASH MANAGER - INTEGRATED FROM STATIC VERSION
====================================================================================================================
*/

// Spatial hash manager for efficient neighbor queries (from static version)
struct spatial_hash_manager {
    float theta_cell_size = 2.0f;  // Will be set based on coarse grid parameters
    float xs_cell_size = 10.0f;    // Will be set based on coarse grid parameters
    
    // Initialize cell sizes based on coarse grid spacing
    void initialize_cell_sizes(float theta_step, float xs_step) {
        theta_cell_size = theta_step / 5.0f;  // Smaller cells for better granularity
        xs_cell_size = xs_step / 3.5f;        // Smaller cells for better granularity
    }
    
    // Hash table: (gg, cell_theta, cell_xs) -> set of job_keys in that cell (O(1) removal)
    std::unordered_map<std::tuple<int, int, int>, std::unordered_set<std::tuple<int, float, float>>> grid;
    
    std::tuple<int, int, int> get_cell_key(int gg, float theta, float xs) const {
        int theta_cell = static_cast<int>(std::floor(theta / theta_cell_size));
        int xs_cell = static_cast<int>(std::floor(xs / xs_cell_size));
        return std::make_tuple(gg, theta_cell, xs_cell);
    }
    
    void add_job(const std::tuple<int, float, float>& job_key) {
        auto [gg, theta, xs] = job_key;
        auto cell_key = get_cell_key(gg, theta, xs);
        grid[cell_key].insert(job_key);
    }
    
    void remove_job(const std::tuple<int, float, float>& job_key) {
        auto [gg, theta, xs] = job_key;
        auto cell_key = get_cell_key(gg, theta, xs);
        auto it = grid.find(cell_key);
        if (it != grid.end()) {
            it->second.erase(job_key);
            if (it->second.empty()) {
                grid.erase(it);
            }
        }
    }
    
    std::vector<std::tuple<int, float, float>> find_neighbors(
        const std::tuple<int, float, float>& job_key, 
        float radius_theta, float radius_xs) const {
        
        auto [gg, theta, xs] = job_key;
        std::vector<std::tuple<int, float, float>> neighbors;
        
        // Calculate how many cells to check in each direction
        // Add extra buffer to handle edge cases and floating-point precision
        int cell_radius_theta = static_cast<int>(std::ceil(radius_theta / theta_cell_size)) + 1;
        int cell_radius_xs = static_cast<int>(std::ceil(radius_xs / xs_cell_size)) + 1;
        
        auto [center_gg, center_theta_cell, center_xs_cell] = get_cell_key(gg, theta, xs);
        
        // Check surrounding cells in a rectangular pattern
        for (int dt = -cell_radius_theta; dt <= cell_radius_theta; ++dt) {
            for (int dx = -cell_radius_xs; dx <= cell_radius_xs; ++dx) {
                auto cell_key = std::make_tuple(center_gg, center_theta_cell + dt, center_xs_cell + dx);
                auto cell_it = grid.find(cell_key);
                
                if (cell_it != grid.end()) {
                    for (const auto& neighbor_job_key : cell_it->second) {
                        auto [neighbor_gg, neighbor_theta, neighbor_xs] = neighbor_job_key;
                        
                        // Skip self
                        if (neighbor_job_key == job_key) continue;
                        
                        // Only consider same gg value
                        if (neighbor_gg != gg) continue;
                        
                        // Check if within radius
                        float d_theta = std::abs(neighbor_theta - theta);
                        float d_xs = std::abs(neighbor_xs - xs);
                        
                        if (d_theta <= radius_theta && d_xs <= radius_xs) {
                            neighbors.push_back(neighbor_job_key);
                        }
                    }
                }
            }
        }
        
        return neighbors;
    }
    
    size_t total_jobs() const {
        size_t total = 0;
        for (const auto& [cell, jobs] : grid) {
            total += jobs.size();
        }
        return total;
    }
};

/*
==================================================================================================================
WORKER
====================================================================================================================
*/

struct worker_state {
    actor manager_actor;
    bool use_early_termination;
};

behavior worker_actor(stateful_actor<worker_state>* self, actor manager, bool use_early_termination) {
    self->state().manager_actor = manager;
    self->state().use_early_termination = use_early_termination;
    
    return {
        [=](job_t job) {
            try {
                auto start_time = std::chrono::high_resolution_clock::now();
                double usmin = -1000.0 / std::pow(2.0, job.n);
                double usmax = 0.0;
                double result = 0.0;
                
                try {
                    // Use C++ solver sweep_main function instead of Julia
                    result = sweep_main(
                        job.Ufile,
                        job.Pfile,
                        job.ufile,
                        job.pfile,
                        job.n,
                        usmax,
                        static_cast<double>(job.gg),
                        static_cast<double>(job.xs),
                        static_cast<double>(job.theta),
                        false,
                        self->state().use_early_termination
                    );
                } catch (const std::exception& e) {
                    self->println("WORKER ERROR: C++ solver exception: {}", e.what());
                    result = 1.0;
                } catch (...) {
                    self->println("WORKER ERROR: Unknown C++ solver exception");
                    result = 1.0;
                }
                
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
                double runtime_sec = duration.count() / 1000.0;
                
                anon_mail("result", job.gg, job.theta, job.xs, result, usmin, usmax, runtime_sec, self)
                    .send(self->state().manager_actor);
                    
            } catch (const std::exception& e) {
                self->println("WORKER EXCEPTION: {}", e.what());
            } catch (...) {
                self->println("WORKER UNKNOWN EXCEPTION");
            }
        },
        [=](const std::string& msg) {
            if (msg == "quit") {
                anon_mail("quit", self).send(self->state().manager_actor);
                self->quit();
            }
        }
    };
}

/*
==================================================================================================================
CLIENT
====================================================================================================================
*/

struct remote_state {
    actor manager;
    bool use_early_termination;
};

behavior remote(stateful_actor<remote_state>* self, const std::string& hostname, uint16_t port, bool enable_early_termination) {
    auto server_actor = self->system().middleman().remote_actor(hostname, port);
    if (!server_actor) {
        self->println("Failed to connect to remote actor at {}:{}", hostname, port);
        self->quit();
        return {};
    }
    
    self->println("Connected to remote actor at {}:{}", hostname, port);
    self->state().manager = *server_actor;
    self->state().use_early_termination = enable_early_termination;
    
    // initial handshake: send (self, 0)
    anon_mail(self, 0).send(self->state().manager);
    
    return {
        [=](int actor_number) {
            for (int i = 0; i < actor_number + 1; ++i) {
                auto worker = self->spawn(worker_actor, self->state().manager, self->state().use_early_termination);
                anon_mail(worker).send(self->state().manager);
            }
            self->quit();
        }
    };
}

/*
==================================================================================================================
SERVER - IMPROVED DYNAMIC WITH SPATIAL HASH
====================================================================================================================
*/

// Enhanced dynamic box structure with spatial hash integration
struct enhanced_dynamic_box {
    int gg;
    float theta_min, theta_max;
    float xs_min, xs_max;
    int level;
    bool refinement_checked = false;   // Has this box been checked for refinement?
    bool refinement_complete = false;  // Has this box completed its refinement?
    
    // Get the 4 corner points of this box
    std::vector<std::tuple<float, float>> get_corner_points() const {
        return {
            {theta_min, xs_min}, // Bottom-left
            {theta_max, xs_min}, // Bottom-right
            {theta_min, xs_max}, // Top-left
            {theta_max, xs_max}  // Top-right
        };
    }
    
    // Subdivide box into 4 children
    std::vector<enhanced_dynamic_box> subdivide() const {
        float theta_mid = (theta_min + theta_max) / 2.0f;
        float xs_mid = (xs_min + xs_max) / 2.0f;
        return {
            {gg, theta_min, theta_mid, xs_min, xs_mid, level + 1, false, false},
            {gg, theta_mid, theta_max, xs_min, xs_mid, level + 1, false, false},
            {gg, theta_min, theta_mid, xs_mid, xs_max, level + 1, false, false},
            {gg, theta_mid, theta_max, xs_mid, xs_max, level + 1, false, false}
        };
    }
    
    // Check if box is small enough (for termination criteria)
    bool is_too_small(float min_theta_size, float min_xs_size) const {
        return (theta_max - theta_min) < min_theta_size || 
               (xs_max - xs_min) < min_xs_size;
    }
    
    // Check if point is inside this box
    bool contains_point(float theta, float xs) const {
        return theta >= theta_min && theta <= theta_max && 
               xs >= xs_min && xs <= xs_max;
    }
    
    // Get a unique box ID for hashing/identification
    std::string get_box_id() const {
        std::ostringstream oss;
        oss << "gg" << gg << "_level" << level << "_theta" << std::fixed << std::setprecision(3) 
            << theta_min << "to" << theta_max << "_xs" << xs_min << "to" << xs_max;
        return oss.str();
    }
    
    // NEW: Check if this box needs neighbors for boundary detection
    std::vector<std::tuple<int, float, float>> get_potential_boundary_neighbors(
        const spatial_hash_manager& spatial_mgr) const {
        
        std::vector<std::tuple<int, float, float>> all_neighbors;
        
        // For each corner point, find neighbors within a reasonable radius
        auto corner_points = get_corner_points();
        float search_radius_theta = (theta_max - theta_min) * 1.5f; // 1.5x box size
        float search_radius_xs = (xs_max - xs_min) * 1.5f;
        
        for (auto [theta, xs] : corner_points) {
            auto job_key = std::make_tuple(gg, theta, xs);
            auto neighbors = spatial_mgr.find_neighbors(job_key, search_radius_theta, search_radius_xs);
            
            // Add neighbors that are not corner points of this box
            for (const auto& neighbor : neighbors) {
                auto [n_gg, n_theta, n_xs] = neighbor;
                if (!contains_point(n_theta, n_xs)) {
                    all_neighbors.push_back(neighbor);
                }
            }
        }
        
        // Remove duplicates
        std::sort(all_neighbors.begin(), all_neighbors.end());
        all_neighbors.erase(std::unique(all_neighbors.begin(), all_neighbors.end()), all_neighbors.end());
        
        return all_neighbors;
    }
};

struct enhanced_server_state {
    enum job_status_t { pending, running, done };
    
    // Job tracking: (gg, theta, xs) -> (status, uq_result)
    std::unordered_map<std::tuple<int, float, float>, std::pair<job_status_t, double>> job_results;
    
    // IMPROVED: Enhanced spatial hash manager (from static version)
    spatial_hash_manager spatial_mgr;
    
    // Dynamic box management: box_id -> box
    std::unordered_map<std::string, enhanced_dynamic_box> active_boxes;
    
    std::deque<job_t> pending_jobs;
    bracket_optimizer::uq_cache_t uq_cache;
    std::unordered_set<std::tuple<int, float, float>> scheduled_jobs;
    
    bool use_bracket;
    bool use_early_termination;
    int active_jobs = 0;
    
    // Termination control
    bool termination_initiated = false;
    
    computation_config comp_config;
    system_config sys_config;
    
    // Worker tracking
    std::unordered_set<actor> active_workers;
    std::unordered_set<actor> idle_workers;
    
    // Statistics
    int total_jobs_created = 0;
    int total_jobs_completed = 0;
    int max_refinement_level = 5;  // Maximum refinement depth
    int boxes_refined = 0;         // Number of boxes that underwent refinement
    int refinement_operations = 0; // Number of individual refinement operations
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point last_activity_time;
    int idle_check_counter = 0;  // Counter to track how many times we've been idle
};

void initialize_job_bracket(stateful_actor<enhanced_server_state>* self, job_t& job) {
    if (self->state().use_bracket) {
        bracket_optimizer::initialize_bracket_from_cache(
            job,
            self->state().uq_cache,
            nullptr, // No neighbor provider for dynamic version
            self->state().sys_config.bracket_config
        );
    }
}

// Helper function to determine if bracket updates should be applied
bool should_apply_bracket_update(stateful_actor<enhanced_server_state>* self, const job_t& job) {
    if (!self->state().use_bracket) return false;
    
    // For dynamic version, we can apply bracket updates more liberally
    // since we're not dealing with synchronized levels
    return true;
}

void wake_idle_workers(stateful_actor<enhanced_server_state>* self);

// Create a new job for the given parameters
void create_job(stateful_actor<enhanced_server_state>* self, int gg, float theta, float xs) {
    auto key = std::make_tuple(gg, theta, xs);
    
    // Check if job already exists
    if (self->state().job_results.find(key) != self->state().job_results.end() || 
        self->state().scheduled_jobs.count(key)) {
        return; // Job already exists
    }
    
    // Create the job
    std::string base = "/globalhome/tus210/HPC/quenchin_actor/";
    std::string waved = base + "waves/index_11";
    std::string crit = base + "waves/index_10";
    
    auto Ufile = waved + "/" + std::to_string(gg) + "/U";
    auto Pfile = waved + "/" + std::to_string(gg) + "/p";
    auto ufile = crit + "/" + std::to_string(gg) + "/U";
    auto pfile = crit + "/" + std::to_string(gg) + "/p";
    
    job_t job{gg, theta, xs, 0, Ufile, Pfile, ufile, pfile};
    self->state().pending_jobs.push_back(job);
    self->state().scheduled_jobs.insert(key);
    self->state().job_results[key] = {enhanced_server_state::pending, 0.0};
    self->state().total_jobs_created++;
    
    // IMPROVED: Add to spatial hash for efficient neighbor queries
    self->state().spatial_mgr.add_job(key);
    
    // Update activity time when creating jobs
    self->state().last_activity_time = std::chrono::high_resolution_clock::now();
    
    self->println("Created job: gg={}, theta={:.3f}, xs={:.3f}", gg, theta, xs);
}

// IMPROVED: Enhanced boundary detection using spatial hash
bool detect_immediate_boundary_with_spatial_hash(stateful_actor<enhanced_server_state>* self, 
                                                const enhanced_dynamic_box& box) {
    auto corner_points = box.get_corner_points();
    std::vector<double> computed_uq_values;
    int completed_corners = 0;
    
    // Check corner completion and collect uq values
    for (auto [theta, xs] : corner_points) {
        auto key = std::make_tuple(box.gg, theta, xs);
        auto result_it = self->state().job_results.find(key);
        
        if (result_it != self->state().job_results.end() && 
            result_it->second.first == enhanced_server_state::done) {
            double uq = result_it->second.second;
            computed_uq_values.push_back(uq);
            completed_corners++;
        }
    }
    
    // Need at least 2 completed corners for early detection
    if (completed_corners < 2) {
        return false;
    }
    
    // Check for sign differences - early boundary detection
    bool has_positive = false, has_negative = false;
    for (double uq : computed_uq_values) {
        if (uq >= 0.0) has_positive = true;
        if (uq < 0.0) has_negative = true;
        
        // Early boundary detection with partial data!
        if (has_positive && has_negative) {
            return true;
        }
    }
    
    // If we have all 4 corners but no boundary detected in corners,
    // check nearby neighbors using spatial hash for more sophisticated boundary detection
    if (completed_corners == 4) {
        auto neighbors = box.get_potential_boundary_neighbors(self->state().spatial_mgr);
        
        for (const auto& neighbor_key : neighbors) {
            auto neighbor_result = self->state().job_results.find(neighbor_key);
            if (neighbor_result != self->state().job_results.end() && 
                neighbor_result->second.first == enhanced_server_state::done) {
                
                double neighbor_uq = neighbor_result->second.second;
                
                // Check if any corner has different sign than this neighbor
                for (double corner_uq : computed_uq_values) {
                    if ((corner_uq >= 0.0) != (neighbor_uq >= 0.0)) {
                        self->println("Boundary detected between box corner (uq={:.6f}) and neighbor (uq={:.6f})", 
                                     corner_uq, neighbor_uq);
                        return true;
                    }
                }
            }
        }
    }
    
    return false; // No boundary detected
}

// IMPROVED: Refine a box using spatial hash for efficient neighbor management
void refine_box_with_spatial_hash(stateful_actor<enhanced_server_state>* self, const enhanced_dynamic_box& box) {
    if (box.level >= self->state().max_refinement_level || box.is_too_small(0.1f, 1.0f)) {
        return; // Skip refinement - max level reached or box too small
    }
    
    self->println("Refining box: gg={}, level={}, theta=[{:.3f},{:.3f}], xs=[{:.3f},{:.3f}]", 
                  box.gg, box.level, box.theta_min, box.theta_max, box.xs_min, box.xs_max);
    
    // Create child boxes
    auto children = box.subdivide();
    int new_jobs_created = 0;
    
    for (const auto& child : children) {
        // Add child box to active boxes
        std::string child_id = child.get_box_id();
        self->state().active_boxes[child_id] = child;
        
        // Create jobs for child box corners that don't exist yet
        auto corner_points = child.get_corner_points();
        for (auto [theta, xs] : corner_points) {
            auto key = std::make_tuple(child.gg, theta, xs);
            
            // Create job if it doesn't exist
            if (self->state().job_results.find(key) == self->state().job_results.end() && 
                self->state().scheduled_jobs.count(key) == 0) {
                create_job(self, child.gg, theta, xs);
                new_jobs_created++;
            }
        }
        
        self->println(" Created child box: level={}, theta=[{:.3f},{:.3f}], xs=[{:.3f},{:.3f}]", 
                      child.level, child.theta_min, child.theta_max, child.xs_min, child.xs_max);
    }
    
    // Remove the parent box since it's now refined
    self->state().active_boxes.erase(box.get_box_id());
    self->state().boxes_refined++;
    self->state().refinement_operations++;
    
    // Wake workers for new jobs
    wake_idle_workers(self);
}

// IMPROVED: Dynamic refinement check using spatial hash for efficient box finding
void check_dynamic_refinement_with_spatial_hash(stateful_actor<enhanced_server_state>* self, 
                                               int gg, float theta, float xs) {
    auto point_key = std::make_tuple(gg, theta, xs);
    
    // IMPROVED: Use spatial hash to efficiently find boxes containing this point
    // Instead of maintaining explicit point-to-box mapping, search spatially
    float search_radius_theta = self->state().sys_config.theta_step * 2.0f; // Reasonable search radius
    float search_radius_xs = self->state().sys_config.xs_step * 2.0f;
    
    std::vector<std::string> boxes_to_check;
    
    // Find boxes that could contain this point
    for (auto& [box_id, box] : self->state().active_boxes) {
        if (box.gg == gg && box.contains_point(theta, xs)) {
            boxes_to_check.push_back(box_id);
        }
    }
    
    // Check ALL boxes that contain this point
    for (const std::string& box_id : boxes_to_check) {
        auto active_box_it = self->state().active_boxes.find(box_id);
        if (active_box_it == self->state().active_boxes.end()) {
            continue; // Box no longer active (might have been refined already)
        }
        
        enhanced_dynamic_box& box = active_box_it->second;
        
        // Skip if this box has already been checked for refinement
        if (box.refinement_checked || box.refinement_complete) {
            continue;
        }
        
        self->println("Checking box for boundary: gg={}, level={}, theta=[{:.3f},{:.3f}], xs=[{:.3f},{:.3f}]",
                      box.gg, box.level, box.theta_min, box.theta_max, box.xs_min, box.xs_max);
        
        // IMPROVED: Use enhanced boundary detection with spatial hash
        if (detect_immediate_boundary_with_spatial_hash(self, box)) {
            self->println("*** BOUNDARY DETECTED - Refining box immediately ***");
            box.refinement_checked = true; // Mark as checked to avoid duplicate processing
            refine_box_with_spatial_hash(self, box);
            continue; // Box was refined, move to next box
        } else {
            // Check if all corners are completed for final decision
            auto corner_points = box.get_corner_points();
            int completed_corners = 0;
            
            for (auto [c_theta, c_xs] : corner_points) {
                auto corner_key = std::make_tuple(box.gg, c_theta, c_xs);
                auto result_it = self->state().job_results.find(corner_key);
                
                if (result_it != self->state().job_results.end() && 
                    result_it->second.first == enhanced_server_state::done) {
                    completed_corners++;
                }
            }
            
            if (completed_corners == 4) {
                // All corners completed, no boundary found - mark as complete
                box.refinement_complete = true;
                self->println("No boundary detected - box refinement complete");
            }
        }
    }
}

void wake_idle_workers(stateful_actor<enhanced_server_state>* self) {
    // Don't wake workers if termination has been initiated
    if (self->state().termination_initiated) {
        return;
    }
    
    auto idle_it = self->state().idle_workers.begin();
    while (idle_it != self->state().idle_workers.end() && !self->state().pending_jobs.empty()) {
        auto idle_worker = *idle_it;
        idle_it = self->state().idle_workers.erase(idle_it);
        
        job_t job = self->state().pending_jobs.front();
        self->state().pending_jobs.pop_front();
        ++self->state().active_jobs;
        
        if (should_apply_bracket_update(self, job)) {
            initialize_job_bracket(self, job);
        }
        
        auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
        self->state().job_results[job_key].first = enhanced_server_state::running;
        self->anon_send(idle_worker, std::move(job));
    }
}

// Initialize the coarse 5x5 grid and create initial boxes with spatial hash
void initialize_coarse_grid_and_boxes_with_spatial_hash(stateful_actor<enhanced_server_state>* self) {
    self->println("=== ENHANCED DYNAMIC BOX-BASED REFINEMENT: Initializing coarse grid ===");
    
    // Initialize spatial hash cell sizes based on coarse grid parameters
    self->state().spatial_mgr.initialize_cell_sizes(self->state().sys_config.theta_step, self->state().sys_config.xs_step);
    
    for (int gg : self->state().sys_config.gg_levels) {
        // Create initial 4x4 boxes at level 0 (covering 5x5 grid)
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                float theta_min = -self->state().sys_config.max_theta + 
                                  (2.0f * self->state().sys_config.max_theta * i) / 4.0f;
                float theta_max = -self->state().sys_config.max_theta + 
                                  (2.0f * self->state().sys_config.max_theta * (i+1)) / 4.0f;
                float xs_min = self->state().sys_config.xs_min + 
                               (self->state().sys_config.xs_max - self->state().sys_config.xs_min) * j / 4.0f;
                float xs_max = self->state().sys_config.xs_min + 
                               (self->state().sys_config.xs_max - self->state().sys_config.xs_min) * (j+1) / 4.0f;
                
                enhanced_dynamic_box box{gg, theta_min, theta_max, xs_min, xs_max, 0, false, false};
                std::string box_id = box.get_box_id();
                self->state().active_boxes[box_id] = box;
                
                // Create jobs for all 4 corners of this box
                auto corner_points = box.get_corner_points();
                for (auto [theta, xs] : corner_points) {
                    create_job(self, gg, theta, xs);
                }
            }
        }
    }
    
    self->println("Initialization complete: {} active boxes, {} pending jobs, {} jobs in spatial hash", 
                  self->state().active_boxes.size(), self->state().pending_jobs.size(),
                  self->state().spatial_mgr.total_jobs());
}

behavior enhanced_server(stateful_actor<enhanced_server_state>* self, uint16_t port, bool enable_bracket, bool enable_early_termination) {
    self->state().use_bracket = enable_bracket;
    self->state().use_early_termination = enable_early_termination;
    self->state().start_time = std::chrono::high_resolution_clock::now();
    self->state().last_activity_time = self->state().start_time;
    
    // Initialize system configuration
    self->state().sys_config = system_config::create_static_config(enable_bracket, enable_early_termination);
    self->state().sys_config.finalize();
    
    // Initialize computation configuration
    self->state().comp_config.enable_bracket_optimization = enable_bracket;
    self->state().comp_config.enable_early_termination = enable_early_termination;
    
    if (!self->system().middleman().publish(self, port)) {
        self->println("Failed to publish server on port {}", port);
        self->quit();
    }
    
    self->println("ENHANCED DYNAMIC BOX Server with Spatial Hash on port {}, bracket-enabled={}, early_termination={}", 
                  port, enable_bracket, enable_early_termination);
    
    // TEST MODE: gamma=3 only, max 3 levels
    self->state().sys_config.gg_levels = {3, 2, 1, 0};
    self->println("=== ENHANCED DYNAMIC BOX-BASED REFINEMENT TEST: gamma=3 only, max 3 levels ===");
    
    // Initialize coarse grid and boxes with spatial hash
    initialize_coarse_grid_and_boxes_with_spatial_hash(self);
    
    // Spawn workers (use full 31 workers for production)
    for (int i = 0; i < self->state().sys_config.actor_number; ++i) {
        auto worker = self->spawn(worker_actor, self, self->state().use_early_termination);
        self->state().active_workers.insert(worker);
        
        if (!self->state().pending_jobs.empty()) {
            job_t job = self->state().pending_jobs.front();
            self->state().pending_jobs.pop_front();
            
            if (should_apply_bracket_update(self, job)) {
                initialize_job_bracket(self, job);
            }
            
            auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
            self->state().job_results[job_key].first = enhanced_server_state::running;
            self->anon_send(worker, std::move(job));
            ++self->state().active_jobs;
        } else {
            self->state().idle_workers.insert(worker);
        }
    }     
    return {
        [=](actor remote, int) {
            anon_mail(self->state().sys_config.actor_number).send(remote);
        },
        
        [=](actor worker) {
            self->state().active_workers.insert(worker);
            
            if (!self->state().pending_jobs.empty()) {
                job_t job = self->state().pending_jobs.front();
                self->state().pending_jobs.pop_front();
                
                if (should_apply_bracket_update(self, job)) {
                    initialize_job_bracket(self, job);
                }
                
                auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
                self->state().job_results[job_key].first = enhanced_server_state::running;
                self->anon_send(worker, std::move(job));
                ++self->state().active_jobs;
            } else {
                self->state().idle_workers.insert(worker);
            }
        },
        
        [=](const std::string& msg, int gg, float theta, float xs, double uq, 
            double usmin, double usmax, double runtime_sec, actor worker) {
            try {
                if (msg != "result") return;
                
                // Store result
                auto job_key = std::make_tuple(gg, theta, xs);
                self->state().job_results[job_key] = {enhanced_server_state::done, uq};
                --self->state().active_jobs;
                ++self->state().total_jobs_completed;
                
                // Update activity time
                self->state().last_activity_time = std::chrono::high_resolution_clock::now();
                
                // Update cache if uq < 0 (for bracket optimization)
                if (uq < 0.0) {
                    // Cache key is (gg, theta), value is vector of result_entry(xs, uq)
                    auto cache_key = std::make_tuple(gg, theta);
                    self->state().uq_cache[cache_key].emplace_back(xs, uq);
                }
                
                self->println("Job completed: gg={}, theta={:.3f}, xs={:.3f}, uq={:.6f}, runtime={:.3f}s", 
                              gg, theta, xs, uq, runtime_sec);
                
                // IMPROVED: Dynamic refinement check with spatial hash
                check_dynamic_refinement_with_spatial_hash(self, gg, theta, xs);
                
                // Assign next job to worker or make worker idle
                if (!self->state().pending_jobs.empty()) {
                    job_t next_job = self->state().pending_jobs.front();
                    self->state().pending_jobs.pop_front();
                    
                    if (should_apply_bracket_update(self, next_job)) {
                        initialize_job_bracket(self, next_job);
                    }
                    
                    auto next_job_key = std::make_tuple(next_job.gg, next_job.theta, next_job.xs);
                    self->state().job_results[next_job_key].first = enhanced_server_state::running;
                    self->anon_send(worker, std::move(next_job));
                } else {
                    self->state().idle_workers.insert(worker);
                }
                
                // Reset idle counter when we have activity (job assignment or worker goes idle)
                self->state().idle_check_counter = 0;
                
                // ENHANCED TERMINATION CHECK: Check if computation is truly complete
                if (self->state().pending_jobs.empty() && self->state().active_jobs == 0) {
                    self->state().idle_check_counter++;
                    
                    // Safety timeout: If we've been idle for too long, force termination
                    auto current_time = std::chrono::high_resolution_clock::now();
                    auto idle_duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - self->state().last_activity_time);
                    
                    if (idle_duration.count() > 60 || self->state().idle_check_counter > 10) {  // 60 seconds idle or 10 idle checks
                        self->println("SAFETY TIMEOUT: Forcing termination after {} seconds idle or {} idle checks", 
                                     idle_duration.count(), self->state().idle_check_counter);
                        
                        if (!self->state().termination_initiated) {
                            self->state().termination_initiated = true;
                            
                            auto end_time = std::chrono::high_resolution_clock::now();
                            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - self->state().start_time);
                            
                            self->println("=== ENHANCED DYNAMIC REFINEMENT COMPLETE (TIMEOUT) ===");
                            self->println("Total runtime: {} seconds", duration.count());
                            self->println("Jobs completed: {}/{}", self->state().total_jobs_completed, self->state().total_jobs_created);
                            self->println("Active boxes remaining: {}", self->state().active_boxes.size());
                            
                            // Shutdown all workers
                            for (auto active_worker : self->state().active_workers) {
                                anon_mail("quit").send(active_worker);
                            }
                            for (auto idle_worker : self->state().idle_workers) {
                                anon_mail("quit").send(idle_worker);
                            }
                        }
                        return;
                    }
                    // Check if all boxes have completed their refinement process
                    bool all_boxes_complete = true;
                    int total_boxes = 0;
                    int complete_boxes = 0;
                    int max_level_boxes = 0;
                    
                    for (const auto& [box_id, box] : self->state().active_boxes) {
                        total_boxes++;
                        
                        if (box.refinement_complete) {
                            complete_boxes++;
                        } else if (box.level >= self->state().max_refinement_level) {
                            max_level_boxes++;  // Boxes at max level are effectively complete
                        } else {
                            // Check if all corners are computed for this box
                            auto corner_points = box.get_corner_points();
                            int completed_corners = 0;
                            
                            for (auto [theta, xs] : corner_points) {
                                auto corner_key = std::make_tuple(box.gg, theta, xs);
                                auto result_it = self->state().job_results.find(corner_key);
                                
                                if (result_it != self->state().job_results.end() && 
                                    result_it->second.first == enhanced_server_state::done) {
                                    completed_corners++;
                                }
                            }
                            
                            // If all corners are computed but box hasn't been marked complete,
                            // it might still be refinable
                            if (completed_corners < 4) {
                                all_boxes_complete = false;
                                break;
                            }
                        }
                    }
                    
                    // Consider computation complete if:
                    // 1. All boxes are marked as refinement_complete, OR
                    // 2. All boxes are at max refinement level, OR
                    // 3. We have processed all corners but no more refinement is happening
                    bool computation_complete = (complete_boxes + max_level_boxes >= total_boxes) || all_boxes_complete;
                    
                    self->println("Termination check: {} total boxes, {} complete, {} at max level, computation_complete={}",
                                  total_boxes, complete_boxes, max_level_boxes, computation_complete);
                    
                    if (computation_complete && !self->state().termination_initiated) {
                        self->state().termination_initiated = true;
                        
                        auto end_time = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - self->state().start_time);
                        
                        self->println("=== ENHANCED DYNAMIC REFINEMENT COMPLETE ===");
                        self->println("Total runtime: {} seconds", duration.count());
                        self->println("Jobs completed: {}/{}", self->state().total_jobs_completed, self->state().total_jobs_created);
                        self->println("Boxes refined: {}", self->state().boxes_refined);
                        self->println("Refinement operations: {}", self->state().refinement_operations);
                        self->println("Final spatial hash size: {} jobs", self->state().spatial_mgr.total_jobs());
                        self->println("Active boxes at termination: {}", self->state().active_boxes.size());
                        
                        // Shutdown all workers
                        for (auto active_worker : self->state().active_workers) {
                            anon_mail("quit").send(active_worker);
                        }
                        for (auto idle_worker : self->state().idle_workers) {
                            anon_mail("quit").send(idle_worker);
                        }
                    }
                }
            } catch (const std::exception& e) {
                self->println("SERVER EXCEPTION in result handler: {}", e.what());
            } catch (...) {
                self->println("SERVER UNKNOWN EXCEPTION in result handler");
            }
        },
        
        [=](const std::string& msg, actor worker) {
            if (msg == "quit") {
                self->state().active_workers.erase(worker);
                self->state().idle_workers.erase(worker);
                self->println("Worker quit received. Remaining workers: {} active, {} idle", 
                              self->state().active_workers.size(), self->state().idle_workers.size());
                
                // Check if all workers have quit
                if (self->state().active_workers.empty() && self->state().idle_workers.empty()) {
                    self->println("All workers have quit. Server shutting down.");
                    self->quit();
                }
            }
        }
    };
}

// MAIN
void caf_main(actor_system& system, const config& cfg) {
    scoped_actor self{system};
    
    self->println("Port: {} Host: {}", cfg.port, cfg.host);
    self->println("Server Mode: {}", cfg.server_mode);
    self->println("Early Termination (Sundials): {}", cfg.enable_early_termination);
    
    if (cfg.server_mode) {
        auto server_actor = system.spawn(enhanced_server, cfg.port, cfg.enable_bracket, cfg.enable_early_termination);
        self->println("Enhanced Server actor spawned");
    } else {
        auto remote_actor = system.spawn(remote, cfg.host, cfg.port, cfg.enable_early_termination);
        self->println("Remote actor spawned");
    }
}

CAF_MAIN(io::middleman, caf::id_block::my_project)