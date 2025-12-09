#include "config.h"
#include "system_config.h"
#include "tuple_hash.h" // Hash specializations for tuples
#include "job_structures.h"
#include "bracket_optimizer.h"
#include "box_neighbor_provider.h"
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
                        job.Ufile,  // fast_U_file
                        job.Pfile,  // fast_P_file
                        job.ufile,  // slow_U_file
                        job.pfile,  // slow_P_file
                        job.n,
                        usmax,
                        static_cast<double>(job.gg),
                        static_cast<double>(job.xs),
                        static_cast<double>(job.theta),
                        false,  // verbose = false
                        self->state().use_early_termination  // use early termination flag
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
SERVER
====================================================================================================================
*/

// Dynamic box structure for immediate refinement checking
struct dynamic_refinement_box {
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
    std::vector<dynamic_refinement_box> subdivide() const {
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
};

struct server_state {
    enum job_status_t { pending, running, done };
    
    // Job tracking: (gg, theta, xs) -> (status, uq_result)
    std::unordered_map<std::tuple<int, float, float>, std::pair<job_status_t, double>> job_results;
    
    // Dynamic box management: box_id -> box
    std::unordered_map<std::string, dynamic_refinement_box> active_boxes;
    
    // Point-to-box mapping for fast lookup: (gg, theta, xs) -> set of box_ids
    std::unordered_map<std::tuple<int, float, float>, std::unordered_set<std::string>> point_to_box_map;
    
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
    
    int total_jobs_created = 0;
    int total_jobs_completed = 0;
    int max_refinement_level = 8;  // Maximum refinement depth
    int boxes_refined = 0;         // Number of boxes that underwent refinement
    int refinement_operations = 0; // Number of individual refinement operations
    std::chrono::high_resolution_clock::time_point start_time;
};

void initialize_job_bracket(stateful_actor<server_state>* self, job_t& job) {
    if (self->state().use_bracket) {
        // Create neighbor provider for spatial optimization
        auto neighbor_provider = std::make_shared<bracket_optimizer::BoxNeighborProvider>(
            reinterpret_cast<const std::unordered_map<std::tuple<int, float, float>, std::pair<int, double>>*>(&self->state().job_results),
            &self->state().point_to_box_map,
            reinterpret_cast<const std::unordered_map<std::string, void*>*>(&self->state().active_boxes),
            0.3f,  // theta tolerance - points within 0.3 units are neighbors
            1.5f   // xs tolerance - points within 1.5 units are neighbors
        );
        
        bracket_optimizer::initialize_bracket_from_cache(
            job,
            self->state().uq_cache,
            neighbor_provider, 
            self->state().sys_config.bracket_config
        );
    }
}

// Helper function to determine if bracket updates should be applied
bool should_apply_bracket_update(stateful_actor<server_state>* self, const job_t& job) {
    if (!self->state().use_bracket) return false;
    
    // For dynamic version, we can apply bracket updates more liberally
    // since we're not dealing with synchronized levels
    return true;
}

void wake_idle_workers(stateful_actor<server_state>* self);

// Create a new job for the given parameters
void create_job(stateful_actor<server_state>* self, int gg, float theta, float xs) {
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
    self->state().job_results[key] = {server_state::pending, 0.0};
    self->state().total_jobs_created++;
    
    self->println("Created job: gg={}, theta={:.3f}, xs={:.3f}", gg, theta, xs);
}

// Enhanced early boundary detection - can detect with partial corner data
bool can_detect_boundary_early(stateful_actor<server_state>* self, const dynamic_refinement_box& box) {
    auto corner_points = box.get_corner_points();
    std::vector<double> computed_uq_values;
    int completed_corners = 0;
    
    for (auto [theta, xs] : corner_points) {
        auto key = std::make_tuple(box.gg, theta, xs);
        auto result_it = self->state().job_results.find(key);
        
        if (result_it != self->state().job_results.end() && 
            result_it->second.first == server_state::done) {
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
        
        // Early boundary detection with partial data
        if (has_positive && has_negative) {
            return true;
        }
    }
    
    return false; // No boundary detected with current data
}

// Check if a box has a quenching boundary (full 4-corner check)
bool box_has_boundary_full(stateful_actor<server_state>* self, const dynamic_refinement_box& box) {
    auto corner_points = box.get_corner_points();
    bool has_positive = false, has_negative = false;
    int completed_corners = 0;
    
    for (auto [theta, xs] : corner_points) {
        auto key = std::make_tuple(box.gg, theta, xs);
        auto result_it = self->state().job_results.find(key);
        
        if (result_it != self->state().job_results.end() && 
            result_it->second.first == server_state::done) {
            double uq = result_it->second.second;
            completed_corners++;
            
            if (uq >= 0.0) has_positive = true;
            if (uq < 0.0) has_negative = true;
        }
    }
    
    // Only consider boundary detection if all 4 corners are complete
    return (completed_corners == 4) && has_positive && has_negative;
}

// Check if all corner points of a box are computed (done or scheduled)
bool box_corners_all_computed_or_scheduled(stateful_actor<server_state>* self, const dynamic_refinement_box& box) {
    auto corner_points = box.get_corner_points();
    
    for (auto [theta, xs] : corner_points) {
        auto key = std::make_tuple(box.gg, theta, xs);
        auto result_it = self->state().job_results.find(key);
        
        // Check if job is either completed or running/pending
        if (result_it == self->state().job_results.end() && 
            self->state().scheduled_jobs.count(key) == 0) {
            return false; // Corner not computed and not scheduled
        }
    }
    
    return true;
}

// Check if all corner points of a box are completed (done status)
bool box_corners_all_completed(stateful_actor<server_state>* self, const dynamic_refinement_box& box) {
    auto corner_points = box.get_corner_points();
    
    for (auto [theta, xs] : corner_points) {
        auto key = std::make_tuple(box.gg, theta, xs);
        auto result_it = self->state().job_results.find(key);
        
        if (result_it == self->state().job_results.end() || 
            result_it->second.first != server_state::done) {
            return false;
        }
    }
    
    return true;
}

// Refine a box by creating 4 child boxes and scheduling corner point jobs
void refine_box(stateful_actor<server_state>* self, const dynamic_refinement_box& box) {
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
            
            // Update point-to-box mapping - ADD to set of boxes (don't overwrite!)
            self->state().point_to_box_map[key].insert(child_id);
            
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
    
    // self->println("Box refinement complete: {} new jobs created", new_jobs_created);
    
    // Wake workers for new jobs
    wake_idle_workers(self);
} // Dynamic refinement check - called after each job completion with early detection
void check_dynamic_refinement(stateful_actor<server_state>* self, int gg, float theta, float xs) {
    auto point_key = std::make_tuple(gg, theta, xs);
    
    // Find which box(es) this point belongs to
    auto box_it = self->state().point_to_box_map.find(point_key);
    if (box_it == self->state().point_to_box_map.end()) {
        return; // Point not associated with any box
    }
    
    // Check ALL boxes that contain this point 
    for (const std::string& box_id : box_it->second) {
        auto active_box_it = self->state().active_boxes.find(box_id);
        if (active_box_it == self->state().active_boxes.end()) {
            continue; // Box no longer active 
        }
        
        dynamic_refinement_box& box = active_box_it->second;
        
        // Skip if this box has already been checked for refinement
        if (box.refinement_checked || box.refinement_complete) {
            continue;
        }
        
        // self->println("Checking box for boundary: gg={}, level={}, theta=[{:.3f},{:.3f}], xs=[{:.3f},{:.3f}]",
        //               box.gg, box.level, box.theta_min, box.theta_max, box.xs_min, box.xs_max);
        
        // FIRST: Try early boundary detection with partial corner data
        if (can_detect_boundary_early(self, box)) {
            // self->println("*** EARLY BOUNDARY DETECTED with partial corners - Refining box immediately ***");
            box.refinement_checked = true; // Mark as checked to avoid duplicate processing
            refine_box(self, box);
            continue; // Box was refined, move to next box
        }
        
        // SECOND: Check if all corners are completed for full boundary analysis
        if (!box_corners_all_completed(self, box)) {
            self->println("Box corners not all completed yet, waiting for more results...");
            continue; // Wait for other corners to complete
        }
        
        // THIRD: All corners completed - do final comprehensive boundary check
        box.refinement_checked = true; // Mark as checked to avoid duplicate checks
        if (box_has_boundary_full(self, box)) {
            // self->println("*** FULL BOUNDARY DETECTED with all corners - Refining box ***");
            refine_box(self, box);
        } else {
            // self->println("No boundary detected - box refinement complete");
            box.refinement_complete = true;
        }
    }
} void wake_idle_workers(stateful_actor<server_state>* self) {
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
        self->state().job_results[job_key].first = server_state::running;
        self->anon_send(idle_worker, std::move(job));
    }
} // Initialize the coarse 5x5 grid and create initial boxes
void initialize_coarse_grid_and_boxes(stateful_actor<server_state>* self) {
    self->println("=== DYNAMIC BOX-BASED REFINEMENT: Initializing coarse grid ===");
    
    for (int gg : self->state().sys_config.gg_levels) {
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
                
                dynamic_refinement_box box{gg, theta_min, theta_max, xs_min, xs_max, 0, false, false};
                std::string box_id = box.get_box_id();
                self->state().active_boxes[box_id] = box;
                
                // Create jobs for all 4 corners of this box and update point-to-box mapping
                auto corner_points = box.get_corner_points();
                for (auto [theta, xs] : corner_points) {
                    auto point_key = std::make_tuple(gg, theta, xs);
                    self->state().point_to_box_map[point_key].insert(box_id);
                    create_job(self, gg, theta, xs);
                }
                
                self->println("Created Level 0 box: gg={}, theta=[{:.3f},{:.3f}], xs=[{:.3f},{:.3f}]", 
                              gg, theta_min, theta_max, xs_min, xs_max);
            }
        }
    }
    
    self->println("Initialization complete: {} active boxes, {} pending jobs", 
                  self->state().active_boxes.size(), self->state().pending_jobs.size());
}

behavior server(stateful_actor<server_state>* self, uint16_t port, bool enable_bracket, bool enable_early_termination) {
    self->state().use_bracket = enable_bracket;
    self->state().use_early_termination = enable_early_termination;
    self->state().start_time = std::chrono::high_resolution_clock::now();
    
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
    
    self->println("DYNAMIC BOX Server up on port {}, bracket-enabled={}, early_termination={}", 
                  port, enable_bracket, enable_early_termination);
    
    self->state().sys_config.gg_levels = {3, 2, 1, 0};
    self->println("=== DYNAMIC BOX-BASED REFINEMENT TEST: gamma=3 only, max 3 levels ===");
    
    // Initialize coarse grid and boxes
    initialize_coarse_grid_and_boxes(self);
    
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
            self->state().job_results[job_key].first = server_state::running;
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
                self->state().job_results[job_key].first = server_state::running;
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
                self->state().job_results[job_key] = {server_state::done, uq};
                --self->state().active_jobs;
                ++self->state().total_jobs_completed;
                
                // Update cache if uq < 0 (for bracket optimization)
                if (uq < 0.0) {
                    self->state().uq_cache[{gg, theta}].emplace_back(xs, uq);
                }
                
                self->println("Job completed: gg={}, theta={:.3f}, xs={:.3f}, uq={:.6f}, bracket=[{:.3f},{:.3f}], runtime={:.3f}s", 
                              gg, theta, xs, uq, usmin, usmax, runtime_sec);
                
                check_dynamic_refinement(self, gg, theta, xs);
                
                // Assign next job to worker or make worker idle
                if (!self->state().pending_jobs.empty()) {
                    job_t job = self->state().pending_jobs.front();
                    self->state().pending_jobs.pop_front();
                    ++self->state().active_jobs;
                    
                    if (should_apply_bracket_update(self, job)) {
                        initialize_job_bracket(self, job);
                    }
                    
                    auto next_job_key = std::make_tuple(job.gg, job.theta, job.xs);
                    self->state().job_results[next_job_key].first = server_state::running;
                    self->anon_send(worker, std::move(job));
                } else {
                    self->state().idle_workers.insert(worker);
                }
                
                // TERMINATION CHECK: no pending jobs + no active jobs = terminate
                if (self->state().pending_jobs.empty() && self->state().active_jobs == 0) {
                    self->println("=== SYSTEM TERMINATION: No pending jobs and no active jobs ===");
                    
                    // Print final statistics
                    int total_computed = self->state().job_results.size();
                    int quenching_points = 0;
                    for (const auto& [key, result] : self->state().job_results) {
                        if (result.second < 0.0) quenching_points++;
                    }
                    
                    self->println("Final statistics:");
                    self->println(" Total points computed: {}", total_computed);
                    self->println(" Quenching points: {} ({:.1f}%)", 
                                  quenching_points, 100.0 * quenching_points / total_computed);
                    self->println(" Boxes refined: {}", self->state().boxes_refined);
                    self->println(" Refinement operations: {}", self->state().refinement_operations);
                    self->println(" Active boxes remaining: {}", self->state().active_boxes.size());
                    
                    auto end_time = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::seconds>(
                        end_time - self->state().start_time);
                    self->println("Total runtime: {} seconds", duration.count());
                    
                    // Terminate all workers
                    for (auto worker : self->state().active_workers) {
                        anon_mail("quit").send(worker);
                    }
                    for (auto worker : self->state().idle_workers) {
                        anon_mail("quit").send(worker);
                    }
                    
                    self->state().termination_initiated = true;
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
                    self->println("All workers have quit. Terminating server...");
                    self->quit();
                }
            }
        }
    };
}

/*
==================================================================================================================
MAIN
====================================================================================================================
*/
void caf_main(actor_system& system, const config& cfg) {
    scoped_actor self{system};
    
    self->println("Port: {} Host: {}", cfg.port, cfg.host);
    self->println("Server Mode: {}", cfg.server_mode);
    self->println("Early Termination (Sundials): {}", cfg.enable_early_termination);
    self->println("Bracket Optimization: {}", cfg.enable_bracket);
    
    if (cfg.server_mode) {
        auto server_actor = system.spawn(server, cfg.port, cfg.enable_bracket, cfg.enable_early_termination);
        self->println("Server actor spawned");
    } else {
        auto remote_actor = system.spawn(remote, cfg.host, cfg.port, cfg.enable_early_termination);
        self->println("Remote actor spawned");
    }
}

CAF_MAIN(io::middleman, caf::id_block::my_project)