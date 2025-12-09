#include "config.h" 
#include "system_config.h"
#include "tuple_hash.h"  // Hash specializations for tuples
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
    f.field("gg",     x.gg),
    f.field("theta",  x.theta),
    f.field("xs",     x.xs),
    f.field("n",      x.n),
    f.field("Ufile",  x.Ufile),
    f.field("Pfile",  x.Pfile),
    f.field("ufile",  x.ufile),
    f.field("pfile",  x.pfile)
  );
}

CAF_BEGIN_TYPE_ID_BLOCK(my_project, caf::first_custom_type_id)
CAF_ADD_TYPE_ID(my_project, (job_t))
CAF_END_TYPE_ID_BLOCK(my_project)

/*
==================================================================================================================
WORKER
=====================================================================================================================
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
            job.Ufile,    // fast_U_file
            job.Pfile,    // fast_P_file
            job.ufile,    // slow_U_file
            job.pfile,    // slow_P_file
            job.n,        
            usmax,   
            static_cast<double>(job.gg),    
            static_cast<double>(job.xs),    
            static_cast<double>(job.theta), 
            false,        // verbose = false
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
=====================================================================================================================
*/

struct remote_state {
  actor manager;
  bool use_early_termination;
};

behavior remote(stateful_actor<remote_state>* self,
                const std::string& hostname,
                uint16_t port,
                bool enable_early_termination) {
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
=====================================================================================================================
*/

// Spatial hash manager for efficient neighbor queries
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
    
    std::tuple<int, int, int> get_cell_key(int gg, float theta, float xs) {
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
            it->second.erase(job_key);  // O(1) removal
            if (it->second.empty()) {
                grid.erase(it);  // Remove empty cell
            }
        }
    }
    
    std::vector<std::tuple<int, float, float>> find_neighbors(
        const std::tuple<int, float, float>& job_key, 
        float radius_theta, float radius_xs) {
        
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
                auto check_cell = std::make_tuple(gg, center_theta_cell + dt, center_xs_cell + dx);
                
                auto it = grid.find(check_cell);
                if (it != grid.end()) {
                    for (auto candidate : it->second) {
                        if (candidate == job_key) continue; // skip self
                        
                        auto [cand_gg, cand_theta, cand_xs] = candidate;
                        
                        // Only consider jobs with the same gg value
                        if (cand_gg != gg) continue;
                        
                        // Calculate distances in both dimensions
                        float d_theta = std::abs(cand_theta - theta);
                        float d_xs = std::abs(cand_xs - xs);
                        
                        // Use rectangular radius check 
                        // This naturally creates 8-directional connectivity
                        if (d_theta <= radius_theta && d_xs <= radius_xs) {
                            neighbors.push_back(candidate);
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

// Box-based adaptive mesh refinement structure
struct refinement_box {
    int gg;
    float theta_min, theta_max;
    float xs_min, xs_max;
    int level;
    bool needs_refinement = false;
    bool processing_complete = false;
    
    // Get the 4 corner points of this box (2x2 grid)
    std::vector<std::tuple<float, float>> get_corner_points() const {
        return {
            {theta_min, xs_min},           // Bottom-left
            {theta_max, xs_min},           // Bottom-right
            {theta_min, xs_max},           // Top-left
            {theta_max, xs_max}            // Top-right
        };
    }
    
    // Subdivide box into 4 children
    std::vector<refinement_box> subdivide() const {
        float theta_mid = (theta_min + theta_max) / 2.0f;
        float xs_mid = (xs_min + xs_max) / 2.0f;
        
        return {
            {gg, theta_min, theta_mid, xs_min, xs_mid, level + 1},
            {gg, theta_mid, theta_max, xs_min, xs_mid, level + 1},
            {gg, theta_min, theta_mid, xs_mid, xs_max, level + 1},
            {gg, theta_mid, theta_max, xs_mid, xs_max, level + 1}
        };
    }
    
    // Check if box is small enough (for termination criteria)
    bool is_too_small(float min_theta_size, float min_xs_size) const {
        return (theta_max - theta_min) < min_theta_size || 
               (xs_max - xs_min) < min_xs_size;
    }
};

struct server_state {

  struct job_node {
    job_t job;
    enum status_t { pending, running, done } status;
    double uq; 
    std::unordered_set<std::tuple<int, float, float>> neighbors; 
  };

  // Box-based refinement system
  std::vector<refinement_box> current_level_boxes;
  std::vector<refinement_box> next_level_boxes;
  
  // Cache of computed results: (gg, theta, xs) -> uq
  std::unordered_map<std::tuple<int, float, float>, double> computed_uq;
  
  std::deque<job_t> pending_jobs;
  bracket_optimizer::uq_cache_t uq_cache;
  std::unordered_set<std::tuple<int, float, float>> scheduled_jobs;
  
  // Legacy structures - kept for compatibility but not used in box-based approach
  std::unordered_map<std::tuple<int, float, float>, job_node> job_graph;
  std::unordered_map<std::tuple<int, float, float>, int> refinement_depth;
  spatial_hash_manager spatial_mgr;

  bool use_bracket;
  bool use_early_termination;
  int  active_jobs = 0;
  
  // Termination control
  bool termination_initiated = false;

  computation_config comp_config;
  
  // System configuration - centralized parameters
  system_config sys_config;

  // Worker tracking
  std::unordered_set<actor> active_workers;
  std::unordered_set<actor> idle_workers;
  
  // Job completion tracking
  int total_jobs_created = 0;
  int total_jobs_completed = 0;

  // Box-based refinement control
  int current_refinement_level = 0;               // Current refinement level being processed
  bool level_processing_complete = false;         // Whether current level is complete
  int max_refinement_level = 8;                   // Maximum refinement depth
  
  // Legacy variables - kept for compatibility
  int global_refinement_level = 0;                
  std::unordered_map<int, bool> gg_level_processing;        
  std::unordered_map<int, bool> gg_completely_done;
  
  // Job tracking per gg and refinement level
  std::unordered_map<std::pair<int, int>, int> gg_level_jobs_created;    
  std::unordered_map<std::pair<int, int>, int> gg_level_jobs_completed;  

  // Refinement timing tracking
  int total_refinement_jobs_created = 0;
  double total_refinement_time_ms = 0.0;  // Legacy timing (ms)
  int refinement_operations = 0;
  
  // Enhanced refinement timing (to match dynamic version)
  double total_refinement_time = 0.0;     // Total time spent in refinement operations (seconds)
  int refinement_time_measurements = 0;   // Number of refinement time measurements
  std::vector<double> level_refinement_times; // Time spent on each refinement level
  
  // Comprehensive refinement cycle timing
  std::unordered_map<int, std::chrono::high_resolution_clock::time_point> refinement_level_start_times;
  std::unordered_map<int, double> refinement_level_durations_ms;
  double total_refinement_cycle_time_ms = 0.0;
  int completed_refinement_levels = 0;
  
  std::chrono::high_resolution_clock::time_point start_time;
};


void initialize_job_bracket(stateful_actor<server_state>* self, job_t& job) {
  if (self->state().use_bracket) {
    bracket_optimizer::initialize_bracket_from_cache(
      job,
      self->state().uq_cache,
      nullptr,  // No neighbor provider for static version
      self->state().sys_config.bracket_config
    );
  }
}



// Helper function to determine if bracket updates should be applied
bool should_apply_bracket_update(stateful_actor<server_state>* self, const job_t& job) {
  if (!self->state().use_bracket) return false;
  
  // For static version, apply bracket updates to all jobs like in dynamic version
  // This ensures that bracket optimization happens at dispatch time, using
  // any available cached data from previously completed jobs
  return true;
}

bool is_quenching_boundary(stateful_actor<server_state>* self, const std::tuple<int, float, float>& job_key, const std::tuple<int, float, float>& neighbor_key);
void trigger_boundary_refinement(stateful_actor<server_state>* self, const std::tuple<int, float, float>& job_key, const std::tuple<int, float, float>& neighbor_key);
void schedule_gg_level_jobs(stateful_actor<server_state>* self, int gg);
void check_gg_completion_and_refine(stateful_actor<server_state>* self);
void cleanup_completed_level_neighbors(stateful_actor<server_state>* self, int completed_level);
void create_all_coarse_grid_neighbors(stateful_actor<server_state>* self);
void wake_idle_workers(stateful_actor<server_state>* self);


// Unified neighbor creation function using euclidean distance for all refinement levels (0-7)
void create_unified_neighbors(stateful_actor<server_state>* self, 
                             const std::tuple<int, float, float>& job_key) {
  auto [gg, theta, xs] = job_key;
  int depth = self->state().refinement_depth[job_key];
  
  // Calculate radius based on depth using unified approach
  float effective_theta_spacing = self->state().sys_config.theta_step / std::pow(2.0f, depth * 0.5f);
  float effective_xs_spacing = self->state().sys_config.xs_step / std::pow(2.0f, depth * 0.5f);
  
  int connectivity_range = self->state().sys_config.connectivity_range;
  float base_radius = connectivity_range * std::min(effective_theta_spacing, effective_xs_spacing);
  
  // For level 0, use coarse_grid_radius_multiplier; for deeper levels, use radius_reduction_factor
  float radius_multiplier;
  if (depth == 0) {
    radius_multiplier = self->state().sys_config.coarse_grid_radius_multiplier;
  } else {
    radius_multiplier = std::pow(self->state().sys_config.radius_reduction_factor, depth);
  }
  
  float euclidean_radius = base_radius * std::sqrt(2.0f) * radius_multiplier;
  
  // Use spatial hash to find potential neighbors within radius
  auto around = self->state().spatial_mgr.find_neighbors(job_key, euclidean_radius, euclidean_radius);
  
  // Connect to all jobs within the euclidean radius
  int connections_made = 0;
  for (auto& existing_key : around) {
    auto [existing_gg, existing_theta, existing_xs] = existing_key;
    
    if (existing_gg != gg) continue; // Only same gg
    
    // Calculate euclidean distance in 2D parameter space
    float d_theta = existing_theta - theta;
    float d_xs = existing_xs - xs;
    float euclidean_distance = std::sqrt(d_theta * d_theta + d_xs * d_xs);
    
    // Connect if within circular radius
    if (euclidean_distance <= euclidean_radius) {
      self->state().job_graph[job_key].neighbors.insert(existing_key);
      self->state().job_graph[existing_key].neighbors.insert(job_key);
      connections_made++;
    }
  }
  
  // Debug output for level 0 to show the improved connectivity
  if (depth == 0) {
    // self->println("Level 0 job connected to {} neighbors with radius {:.1f}", connections_made, euclidean_radius);
  }
}

void create_coarse_grid_neighbors(stateful_actor<server_state>* self, 
                                 const std::tuple<int, float, float>& job_key) {
  // Use unified neighbor creation for all levels, including level 0
  create_unified_neighbors(self, job_key);
}

void schedule_job(stateful_actor<server_state>* self, job_t new_job, bool is_refinement = false) {
  
  auto key = std::make_tuple(new_job.gg, new_job.theta, new_job.xs);

  self->state().job_graph[key] = server_state::job_node{new_job, server_state::job_node::pending, 0.0, {}};
  
  self->state().total_jobs_created++;

  self->state().spatial_mgr.add_job(key);

  // Explicitly initialize refinement depth for new jobs
  if (self->state().refinement_depth.find(key) == self->state().refinement_depth.end()) {
      self->state().refinement_depth[key] = 0; // Default to coarse grid job
  }

  // Add to pending_jobs + scheduled_jobs (ready for workers + prevents duplicates)
  if (is_refinement) {
    self->state().pending_jobs.push_back(new_job);
  } else {
    self->state().pending_jobs.push_back(new_job);
  }
  self->state().scheduled_jobs.insert(key);
  
  int depth = self->state().refinement_depth[key];
}

void check_gg_completion_and_refine(stateful_actor<server_state>* self) {
  
  bool all_gg_completed_current_level = true;
  int current_level = self->state().global_refinement_level;
  
  for (int gg : self->state().sys_config.gg_levels) {
    // Skip if completely done
    if (self->state().gg_completely_done[gg]) continue;
    
    auto level_key = std::make_pair(gg, current_level);
    int level_created = self->state().gg_level_jobs_created[level_key];
    int level_completed = self->state().gg_level_jobs_completed[level_key];
    
    if (level_completed < level_created || level_created == 0) {
      all_gg_completed_current_level = false;
    }
  }
  
  // If ALL gg levels completed current level, advance ALL to next level
  if (all_gg_completed_current_level) {
    // Check if we can advance to next refinement level
    int next_level = current_level + 1;
    if (next_level > self->state().sys_config.max_refinement_level) {
      // All levels reached maximum refinement
      
      // Mark all as completely done
      for (int gg : self->state().sys_config.gg_levels) {
        self->state().gg_completely_done[gg] = true;
      }
    } else {
      // Record end time for previous refinement level
      if (current_level > 0) {
        auto current_time = std::chrono::high_resolution_clock::now();
        if (self->state().refinement_level_start_times.find(current_level) != self->state().refinement_level_start_times.end()) {
          auto level_duration = std::chrono::duration_cast<std::chrono::microseconds>(
            current_time - self->state().refinement_level_start_times[current_level]);
          double level_time_ms = level_duration.count() / 1000.0;
          self->state().refinement_level_durations_ms[current_level] = level_time_ms;
          self->state().total_refinement_cycle_time_ms += level_time_ms;
          self->state().completed_refinement_levels++;
        }
      }
      
      // Advance global refinement level
      self->state().global_refinement_level = next_level;
      
      // Start timing the new refinement level
      self->state().refinement_level_start_times[next_level] = std::chrono::high_resolution_clock::now();
      
      // Create refinement jobs for ALL gg levels at the new level
      int total_refinement_jobs = 0;
      
      for (int gg : self->state().sys_config.gg_levels) {
        if (!self->state().gg_completely_done[gg]) {
          // Mark as processing next level
          self->state().gg_level_processing[gg] = true;
          
          // Find boundary pairs for this gg at the completed level
          std::set<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> boundary_pairs;
          std::set<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> processed_pairs;
          
          // Scan all completed jobs in this gg level for boundaries
          int completed_jobs_in_gg = 0;
          int total_neighbors = 0;
          for (auto& [job_key, job_node] : self->state().job_graph) {
            auto [job_gg, job_theta, job_xs] = job_key;
            if (job_gg != gg || job_node.status != server_state::job_node::done) continue;
            
            completed_jobs_in_gg++;
            total_neighbors += job_node.neighbors.size();
            
            // Create a copy of neighbors to avoid iterator invalidation
            auto neighbors_copy = job_node.neighbors;
            
            // Check each neighbor for boundary conditions
            for (auto neighbor_key : neighbors_copy) {
              auto& neighbor_node = self->state().job_graph[neighbor_key];
              auto [neighbor_gg, neighbor_theta, neighbor_xs] = neighbor_key;
              
              // Only consider neighbors in same gg level that are completed
              if (neighbor_gg == gg && neighbor_node.status == server_state::job_node::done) {
                
                // Pre-filter boundary pairs that have already reached max depth
                int depth1 = self->state().refinement_depth[job_key];
                int depth2 = self->state().refinement_depth[neighbor_key];
                int max_parent_depth = std::max(depth1, depth2);
                
                if (max_parent_depth >= self->state().sys_config.max_refinement_level) {
                  continue; // Skip - this pair has already reached maximum refinement depth
                }
                
                if (is_quenching_boundary(self, job_key, neighbor_key)) {
                  // Create ordered pair to avoid duplicates
                  auto pair1 = std::make_pair(job_key, neighbor_key);
                  auto pair2 = std::make_pair(neighbor_key, job_key);
                  
                  // Only add if we haven't seen this pair in either direction
                  if (boundary_pairs.count(pair1) == 0 && boundary_pairs.count(pair2) == 0) {
                    boundary_pairs.insert(pair1);
                  }
                }
              }
            }
          }
          
        //   self->println("Found {} boundary pairs for gg {} level {} → {}", 
        //                 boundary_pairs.size(), gg, current_level, next_level);
          
          // Create refinement jobs for all boundary pairs
          int gg_refinement_jobs = 0;
          auto next_level_key = std::make_pair(gg, next_level);
          
          for (auto& boundary_pair : boundary_pairs) {
            auto job_key = boundary_pair.first;
            auto neighbor_key = boundary_pair.second;
            
            // Ensure we don't process the same pair twice
            if (processed_pairs.count(boundary_pair) == 0) {
              auto reverse_pair = std::make_pair(neighbor_key, job_key);
              processed_pairs.insert(boundary_pair);
              processed_pairs.insert(reverse_pair);  // Mark both directions as processed
              
              // Count jobs before refinement
              int jobs_before = self->state().total_jobs_created;
              
              trigger_boundary_refinement(self, job_key, neighbor_key);
              
              // Count jobs after refinement
              if (self->state().total_jobs_created > jobs_before) {
                gg_refinement_jobs++;
              }
            }
          }
          
          // Update job counts for this gg at next level
          self->state().gg_level_jobs_created[next_level_key] = gg_refinement_jobs;
          self->state().gg_level_jobs_completed[next_level_key] = 0;
          total_refinement_jobs += gg_refinement_jobs;
          
          self->println("=== GG {} advanced to level {} with {} jobs created ===", 
                        gg, next_level, gg_refinement_jobs);
          
          // If no jobs created for this gg, mark as completely done
          if (gg_refinement_jobs == 0) {
            self->state().gg_completely_done[gg] = true;
            // self->println("=== GG {} COMPLETELY FINISHED - no boundary refinement possible ===", gg);
          }
          
          self->state().gg_level_processing[gg] = false;
        }
      }
      
      if (total_refinement_jobs > 0) {
        self->println("=== SYNCHRONIZED REFINEMENT: Advanced ALL gg levels to level {} with {} total jobs ===", 
                      next_level, total_refinement_jobs);
        
        // Wake workers for new jobs
        wake_idle_workers(self);
      } else {
        // No refinement jobs created for any gg - all are done
        self->println("=== ALL GG LEVELS COMPLETELY FINISHED - no more refinement possible ===");
        for (int gg : self->state().sys_config.gg_levels) {
          self->state().gg_completely_done[gg] = true;
        }
      }
      
      // Clean up obsolete neighbor relationships from the completed level
      // (do this AFTER boundary detection and refinement job creation)
      cleanup_completed_level_neighbors(self, current_level);
    }
  }
  
  // Check if ALL gg levels are done
  bool all_done = true;
  for (int gg : self->state().sys_config.gg_levels) {
    if (!self->state().gg_completely_done[gg]) {
      all_done = false;
      break;
    }
  }
  
  if (all_done) {
    // self->println("=== ALL GG LEVELS COMPLETED ===");
    
    // Shutdown all workers
    for (auto worker : self->state().active_workers) {
      anon_mail("quit").send(worker);
    }
    for (auto worker : self->state().idle_workers) {
      anon_mail("quit").send(worker);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - self->state().start_time);
    self->println("Total runtime: {} seconds", duration.count());
    self->println("Jobs completed: {}/{}", self->state().total_jobs_completed, self->state().total_jobs_created);
    
    // Refinement timing summary
    self->println("=== REFINEMENT TIMING SUMMARY (STATIC) ===");
    
    // Job creation overhead timing (original measurement)
    self->println("Refinement Job Creation Overhead:");
    self->println("  Total refinement operations: {}", self->state().refinement_operations);
    self->println("  Total refinement jobs created: {}", self->state().total_refinement_jobs_created);
    self->println("  Job creation time: {:.3f} ms", self->state().total_refinement_time_ms);
    if (self->state().refinement_operations > 0) {
        double avg_refinement_time = self->state().total_refinement_time_ms / self->state().refinement_operations;
        self->println("  Average job creation time: {:.3f} ms per operation", avg_refinement_time);
    }
    
    // Enhanced refinement timing (unified with dynamic version)
    self->println("Enhanced Refinement Timing:");
    self->println("  Total refinement time: {:.6f}s ({} level operations)", 
                 self->state().total_refinement_time, self->state().refinement_time_measurements);
    
    if (self->state().refinement_time_measurements > 0) {
        double avg_refinement_time = self->state().total_refinement_time / self->state().refinement_time_measurements;
        self->println("  Average refinement time per level: {:.6f}s", avg_refinement_time);
    }
    
    // Per-level timing breakdown
    for (size_t i = 0; i < self->state().level_refinement_times.size(); ++i) {
        self->println("  Level {} refinement time: {:.6f}s", i, self->state().level_refinement_times[i]);
    }
    
    // Percentage calculations
    double job_creation_percentage = (duration.count() > 0) ? 
        (self->state().total_refinement_time_ms / 1000.0) / duration.count() * 100.0 : 0.0;
    double enhanced_refinement_percentage = (duration.count() > 0) ? 
        self->state().total_refinement_time / duration.count() * 100.0 : 0.0;
        
    self->println("Timing as Percentage of Total Runtime:");
    self->println("  Job creation overhead: {:.3f}%", job_creation_percentage);
    self->println("  Enhanced refinement overhead: {:.2f}%", enhanced_refinement_percentage);
    
    self->state().termination_initiated = true;
    self->quit();
  } else {
    wake_idle_workers(self);
  }
}

void wake_idle_workers(stateful_actor<server_state>* self) {
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
    
    // Apply bracket optimization for jobs at dispatch time (like dynamic version)
    if (should_apply_bracket_update(self, job)) {
      initialize_job_bracket(self, job);
    }
    
    auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
    self->state().job_graph[job_key].status = server_state::job_node::running;
    self->anon_send(idle_worker, std::move(job));
  }
}

// Schedule jobs for a specific gg level
void schedule_gg_level_jobs(stateful_actor<server_state>* self, int gg) {
  float max_theta = self->state().sys_config.max_theta;
  int theta_step = self->state().sys_config.theta_step;
  const float xs_min = self->state().sys_config.xs_min;
  const float xs_max = self->state().sys_config.xs_max;
  int NDNS = self->state().sys_config.NDNS;
  
  std::vector<float> xs_values;
  for (int i = 0; i < NDNS; ++i) {
    float xs = xs_min + i * self->state().sys_config.xs_step;  // Use pre-calculated xs_step
    xs_values.push_back(xs);
  }

  std::vector<int> all_theta;
  for (float theta = -max_theta; theta <= max_theta; theta += theta_step)
    all_theta.push_back(theta);

  std::string base = "/globalhome/tus210/HPC/quenchin_actor/";
  std::string waved = base + "waves/index_11";
  std::string crit = base + "waves/index_10";
  auto Ufile = waved + "/" + std::to_string(gg) + "/U";
  auto Pfile = waved + "/" + std::to_string(gg) + "/p";
  auto ufile = crit + "/" + std::to_string(gg) + "/U";
  auto pfile = crit + "/" + std::to_string(gg) + "/p";

  int jobs_for_gg = 0;
  for (float xs : xs_values) {
    for (float theta : all_theta) {
      auto job_key = std::make_tuple(gg, theta, xs);
      if (self->state().scheduled_jobs.count(job_key)) continue;
      
      self->state().scheduled_jobs.insert(job_key);

      job_t job{
        gg, theta, xs, 0,  // n=0 gives usmin=-1000, usmax=0
        Ufile, Pfile, ufile, pfile
      };
      schedule_job(self, job);
      jobs_for_gg++;
    }
  }
  
  // Initialize tracking for this gg level (level 0) using global refinement level
  auto level_0_key = std::make_pair(gg, 0);
  self->state().gg_level_jobs_created[level_0_key] = jobs_for_gg;
  self->state().gg_level_jobs_completed[level_0_key] = 0;
  self->state().gg_completely_done[gg] = false;
}

void bracket_initialize_from_cache(job_t& job,
  const bracket_optimizer::uq_cache_t& uq_cache,
  float L)
{
  // Empty implementation - functionality moved to bracket_optimizer module
  // This function is kept for compatibility but actual implementation
  // is handled by initialize_job_bracket which calls bracket_optimizer::initialize_bracket_from_cache
}

bool is_quenching_boundary(stateful_actor<server_state>* self,
                           const std::tuple<int, float, float>& job_key,
                           const std::tuple<int, float, float>& neighbor_key) {
  // Detects boundaries between two completed jobs for adaptive refinement
  
  auto& job_node = self->state().job_graph[job_key];
  auto& neighbor_node = self->state().job_graph[neighbor_key];
  
  double uq1 = job_node.uq;
  double uq2 = neighbor_node.uq;
  
  // Primary boundary criterion: one quenches, one doesn't
  bool job_quenches = (uq1 < 0.0);
  bool neighbor_quenches = (uq2 < 0.0);
  bool is_boundary = false;
  
  if (job_quenches != neighbor_quenches) {
    // Debug output for true boundaries
    // self->println("*** TRUE BOUNDARY DETECTED: ({},{},{}) uq={} vs ({},{},{}) uq={} ***", 
    //               std::get<0>(job_key), std::get<1>(job_key), std::get<2>(job_key), uq1,
    //               std::get<0>(neighbor_key), std::get<1>(neighbor_key), std::get<2>(neighbor_key), uq2);
    is_boundary = true;
  }
  
  // OPTIMIZATION: Remove this neighbor pair from both jobs' neighbor lists
  // since we've now processed them and won't need to check again
  // If boundary found: the refinement job created between them will serve as the bridge
  // If no boundary: no need to check this pair again
  self->state().job_graph[job_key].neighbors.erase(neighbor_key);
  self->state().job_graph[neighbor_key].neighbors.erase(job_key);
  
  return is_boundary;
}

void trigger_boundary_refinement(stateful_actor<server_state>* self,
                                const std::tuple<int, float, float>& job_key,
                                const std::tuple<int, float, float>& neighbor_key) {
  // Start timing this refinement operation
  auto refinement_start_time = std::chrono::high_resolution_clock::now();
  
  int max_depth = self->state().sys_config.max_refinement_level;  // Use dynamic max refinement level
  auto [gg1, theta1, xs1] = job_key;
  auto [gg2, theta2, xs2] = neighbor_key;

  // Only refine within the same gg value
  if (gg1 != gg2) {
    return;  // Skip silently - different gg values
  }

  // Check refinement depth for both jobs
  int depth1 = self->state().refinement_depth[job_key];
  int depth2 = self->state().refinement_depth[neighbor_key];
  
  // Check if we would create jobs beyond max_depth
  int max_parent_depth = std::max(depth1, depth2);
  if (max_parent_depth >= max_depth) {
    // Skip silently - max depth reached (this is expected and normal)
    return;
  }

  // Calculate midpoint between the two boundary jobs
  float center_theta = (theta1 + theta2) / 2.0f;
  float center_xs = (xs1 + xs2) / 2.0f;

  // Check bounds - use dynamic max_theta instead of hard-coded values
  if (center_theta < -self->state().sys_config.max_theta || center_theta > self->state().sys_config.max_theta ||
      center_xs < self->state().sys_config.xs_min || center_xs > self->state().sys_config.xs_max) {
    return;  // Skip silently - midpoint out of bounds
  }
  
  auto new_job_key = std::make_tuple(gg1, center_theta, center_xs);
  if (self->state().job_graph.count(new_job_key) > 0) {
    return;  // Skip silently - midpoint job already exists
  }

  // Calculate the refinement depth for new job (increment from max parent depth)
  int new_depth = max_parent_depth + 1;
  
  // DISTANCE-BASED THRESHOLD: Prevent over-refinement by checking minimum distances
  // Scale-based thresholds that decrease with refinement depth
  // COMMENTED OUT: Allow nodes to get closer for better refinement resolution
  /*
  float min_theta_distance = self->state().theta_step / std::pow(2.0f, new_depth + 1);  // Exponential decrease
  float min_xs_distance = self->state().xs_step / std::pow(2.0f, new_depth + 1);        // Exponential decrease
  
  // Ensure minimum absolute thresholds for numerical stability
  min_theta_distance = std::max(min_theta_distance, 0.5f);   // At least 0.5 degrees
  min_xs_distance = std::max(min_xs_distance, 1.0f);         // At least 1.0 xs units
  
  // Check if potential refinement point is too close to existing jobs
  // Use spatial hash to efficiently find nearby jobs
  auto nearby_jobs = self->state().spatial_mgr.find_neighbors(new_job_key, 
                                                              min_theta_distance * 2.0f, 
                                                              min_xs_distance * 2.0f);
  
  for (auto& existing_key : nearby_jobs) {
    auto [existing_gg, existing_theta, existing_xs] = existing_key;
    
    // Calculate distances to existing job
    float d_theta = std::abs(existing_theta - center_theta);
    float d_xs = std::abs(existing_xs - center_xs);
    
    // Skip refinement if too close to an existing job
    if (d_theta < min_theta_distance && d_xs < min_xs_distance) {
      // self->println("REFINEMENT: Skipping - too close to existing job ({},{:.3f},{:.3f}). Distance=({:.3f},{:.3f}) < threshold=({:.3f},{:.3f})", 
      //               existing_gg, existing_theta, existing_xs, d_theta, d_xs, min_theta_distance, min_xs_distance);
      return;  // Skip silently - too close to existing job
    }
  }
  */

  // Create new job
  std::string base = "/globalhome/tus210/HPC/quenchin_actor/";
  std::string waved = base + "waves/index_11";
  std::string crit = base + "waves/index_10";

  job_t new_job{
      gg1, center_theta, center_xs, 0,  // n=0 gives usmin=-1000, usmax=0
      waved + "/" + std::to_string(gg1) + "/U",
      waved + "/" + std::to_string(gg1) + "/p",
      crit + "/" + std::to_string(gg1) + "/U",
      crit + "/" + std::to_string(gg1) + "/p"
  };
  
  self->state().refinement_depth[new_job_key] = new_depth;
  
  // Bracket optimization is now handled at dispatch time for consistency
  
  // Schedule the job to add it to spatial hash and job graph (add to END of queue)
  schedule_job(self, new_job, true);  // true = is_refinement

  // Track refinement jobs per gg level (use global refinement level)
  int current_level = self->state().global_refinement_level;
  auto level_key = std::make_pair(gg1, current_level);
  self->state().gg_level_jobs_created[level_key]++;

  // Create neighbors using unified approach - ensures consistent neighbor creation for all levels
  int initial_neighbors = self->state().job_graph[new_job_key].neighbors.size();
  create_unified_neighbors(self, new_job_key);
  int connections_made = self->state().job_graph[new_job_key].neighbors.size() - initial_neighbors;
  
  // self->println("Refinement job connected to {} neighbors within euclidean radius {:.1f}", 
  //              connections_made, euclidean_radius);
  
  // Record timing for this refinement operation
  auto refinement_end_time = std::chrono::high_resolution_clock::now();
  auto refinement_duration = std::chrono::duration_cast<std::chrono::microseconds>(refinement_end_time - refinement_start_time);
  double refinement_time_ms = refinement_duration.count() / 1000.0;
  
  self->state().total_refinement_time_ms += refinement_time_ms;
  self->state().total_refinement_jobs_created++;
  self->state().refinement_operations++;
  
//   self->println("REFINEMENT: Created job (gg={} global_level {} job #{}): gg={}, theta={}, xs={} (depth={}) with {} neighbors in {:.3f}ms", 
//                 gg1, current_level, self->state().gg_level_jobs_created[level_key], gg1, center_theta, center_xs, new_depth, connections_made, refinement_time_ms);
}

// Function to create complete coarse grid neighbor map after all jobs are registered
void create_all_coarse_grid_neighbors(stateful_actor<server_state>* self) {
  // Legacy function - kept for compatibility but not used in box-based approach
  self->println("Box-based approach: skipping legacy neighbor creation");
}

// ============================================================================
// BOX-BASED REFINEMENT FUNCTIONS
// ============================================================================

// Check if a box contains quenching boundary (using its 4 corner points)
bool box_has_boundary(stateful_actor<server_state>* self, const refinement_box& box) {
    auto corner_points = box.get_corner_points();
    
    bool has_positive = false, has_negative = false;
    int computed_corners = 0;
    
    for (auto [theta, xs] : corner_points) {
        auto key = std::make_tuple(box.gg, theta, xs);
        
        if (self->state().computed_uq.count(key)) {
            double uq = self->state().computed_uq[key];
            computed_corners++;
            
            if (uq >= 0.0) has_positive = true;
            if (uq < 0.0) has_negative = true;
            
            // Early exit if boundary found
            if (has_positive && has_negative) {
                return true;
            }
        }
    }
    
    // Only check boundary if we have data for all 4 corners
    return (computed_corners == 4) && has_positive && has_negative;
}

// Create initial 5x5 grid and organize into boxes
void initialize_coarse_grid_and_boxes(stateful_actor<server_state>* self) {
    self->state().current_level_boxes.clear();
    
    for (int gg : self->state().sys_config.gg_levels) {
        // Create 5x5 initial grid jobs
        std::vector<std::tuple<float, float>> grid_points;
        
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                float theta = -self->state().sys_config.max_theta + 
                             (2.0f * self->state().sys_config.max_theta * i) / 4.0f;
                float xs = self->state().sys_config.xs_min + 
                          (self->state().sys_config.xs_max - self->state().sys_config.xs_min) * j / 4.0f;
                grid_points.emplace_back(theta, xs);
            }
        }
        
        // Create jobs for all 25 grid points
        for (auto [theta, xs] : grid_points) {
            auto key = std::make_tuple(gg, theta, xs);
            
            if (self->state().scheduled_jobs.find(key) == self->state().scheduled_jobs.end()) {
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
                
                self->println("  Created initial grid job: gamma={}, theta={:.2f}, xs={:.2f}", 
                             gg, theta, xs);
            }
        }
        
        // Create 4x4 = 16 boxes from the 5x5 grid
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
                
                refinement_box box{gg, theta_min, theta_max, xs_min, xs_max, 0, false, false};
                self->state().current_level_boxes.push_back(box);
                
                self->println("  Created Level 0 box: gamma={}, theta=[{:.1f},{:.1f}], xs=[{:.1f},{:.1f}]", 
                             gg, theta_min, theta_max, xs_min, xs_max);
            }
        }
    }
    
    self->println("Initial grid: 25 jobs created, {} boxes at level 0", 
                  self->state().current_level_boxes.size());
}

// Process jobs for boxes that need new corner points
void process_current_level(stateful_actor<server_state>* self) {
    self->println("=== Processing refinement level {} with {} boxes ===", 
                  self->state().current_refinement_level,
                  self->state().current_level_boxes.size());
    
    int jobs_created = 0;
    
    // For level > 0, we already have the boxes and only need to create jobs for missing corner points
    if (self->state().current_refinement_level > 0) {
        for (const auto& box : self->state().current_level_boxes) {
            if (box.level == self->state().current_refinement_level) {
                auto corner_points = box.get_corner_points();
                
                for (auto [theta, xs] : corner_points) {
                    auto key = std::make_tuple(box.gg, theta, xs);
                    
                    // Only create job if not already computed or scheduled
                    if (self->state().computed_uq.find(key) == self->state().computed_uq.end() &&
                        self->state().scheduled_jobs.find(key) == self->state().scheduled_jobs.end()) {
                        
                        std::string base = "/globalhome/tus210/HPC/quenchin_actor/";
                        std::string waved = base + "waves/index_11";
                        std::string crit = base + "waves/index_10";
                        auto Ufile = waved + "/" + std::to_string(box.gg) + "/U";
                        auto Pfile = waved + "/" + std::to_string(box.gg) + "/p";
                        auto ufile = crit + "/" + std::to_string(box.gg) + "/U";
                        auto pfile = crit + "/" + std::to_string(box.gg) + "/p";
                        
                        job_t job{box.gg, theta, xs, 0, Ufile, Pfile, ufile, pfile};
                        
                        self->println("  Created job: gamma={}, theta={:.2f}, xs={:.2f} (level {})", 
                                     box.gg, theta, xs, box.level);
                        
                        self->state().pending_jobs.push_back(job);
                        self->state().scheduled_jobs.insert(key);
                        jobs_created++;
                    }
                }
            }
        }
    }
    
    self->println("Level {}: Created {} new jobs for corner points, total pending: {}", 
                  self->state().current_refinement_level, jobs_created, self->state().pending_jobs.size());
    
    // Wake idle workers to start processing
    wake_idle_workers(self);
}

// Check if current level is complete and advance to next level
void check_level_completion_and_advance(stateful_actor<server_state>* self) {
    // Don't process completion checks if termination has been initiated
    if (self->state().termination_initiated) {
        return;
    }
    
    // Check if all corner points for current level boxes are computed
    int corners_needed = 0, corners_completed = 0;
    
    for (const auto& box : self->state().current_level_boxes) {
        if (box.level == self->state().current_refinement_level) {
            auto corner_points = box.get_corner_points();
            
            for (auto [theta, xs] : corner_points) {
                auto key = std::make_tuple(box.gg, theta, xs);
                corners_needed++;
                
                if (self->state().computed_uq.count(key)) {
                    corners_completed++;
                }
            }
        }
    }
    
    if (corners_completed < corners_needed) {
        return; // Level not complete yet
    }
    
    self->println("=== Level {} COMPLETE: {}/{} corner points computed ===", 
                  self->state().current_refinement_level, corners_completed, corners_needed);
    
    // Start timing the complete refinement cycle for this level (boundary detection + box creation)
    auto level_refinement_start = std::chrono::high_resolution_clock::now();
    
    // Detect boundaries and create next level boxes
    std::vector<refinement_box> next_level_boxes;
    int boxes_with_boundaries = 0;
    int new_jobs_created = 0;
    
    for (const auto& box : self->state().current_level_boxes) {
        if (box.level == self->state().current_refinement_level) {
            
            if (box_has_boundary(self, box) && 
                box.level < self->state().max_refinement_level &&
                !box.is_too_small(0.1f, 1.0f)) {
                
                // Subdivide this box and create jobs for new corner points
                auto children = box.subdivide();
                
                for (const auto& child_box : children) {
                    next_level_boxes.push_back(child_box);
                    
                    // Create jobs for any new corner points not already computed
                    auto corner_points = child_box.get_corner_points();
                    for (auto [theta, xs] : corner_points) {
                        auto key = std::make_tuple(child_box.gg, theta, xs);
                        
                        if (self->state().computed_uq.find(key) == self->state().computed_uq.end() &&
                            self->state().scheduled_jobs.find(key) == self->state().scheduled_jobs.end()) {
                            
                            std::string base = "/globalhome/tus210/HPC/quenchin_actor/";
                            std::string waved = base + "waves/index_11";
                            std::string crit = base + "waves/index_10";
                            auto Ufile = waved + "/" + std::to_string(child_box.gg) + "/U";
                            auto Pfile = waved + "/" + std::to_string(child_box.gg) + "/p";
                            auto ufile = crit + "/" + std::to_string(child_box.gg) + "/U";
                            auto pfile = crit + "/" + std::to_string(child_box.gg) + "/p";
                            
                            job_t job{child_box.gg, theta, xs, 0, Ufile, Pfile, ufile, pfile};
                            
                            self->state().pending_jobs.push_back(job);
                            self->state().scheduled_jobs.insert(key);
                            new_jobs_created++;
                        }
                    }
                }
                
                boxes_with_boundaries++;
                
                self->println("Box subdivided: gamma={}, level {} -> {}, theta=[{:.2f},{:.2f}], xs=[{:.2f},{:.2f}]",
                             box.gg, box.level, box.level + 1,
                             box.theta_min, box.theta_max, box.xs_min, box.xs_max);
            }
        } else if (box.level > self->state().current_refinement_level) {
            // Keep boxes from higher levels
            next_level_boxes.push_back(box);
        }
    }
    
    // Measure refinement time for this level
    auto level_refinement_end = std::chrono::high_resolution_clock::now();
    auto level_duration = std::chrono::duration_cast<std::chrono::microseconds>(level_refinement_end - level_refinement_start);
    double level_time_sec = level_duration.count() / 1e6;
    
    // Update timing statistics
    self->state().total_refinement_time += level_time_sec;
    self->state().refinement_time_measurements++;
    self->state().level_refinement_times.push_back(level_time_sec);
    
    // Update current level
    self->state().current_level_boxes = std::move(next_level_boxes);
    self->state().current_refinement_level++;
    
    self->println("Level {}: {} boxes had boundaries, {} new boxes created, {} new corner jobs, refinement time: {:.6f}s",
                  self->state().current_refinement_level - 1,
                  boxes_with_boundaries,
                  self->state().current_level_boxes.size(),
                  new_jobs_created,
                  level_time_sec);
    
    if (self->state().current_level_boxes.empty()) {
        self->println("=== REFINEMENT COMPLETE - No more boundaries found ===");
        
        int total_computed = self->state().computed_uq.size();
        int quenching_points = 0;
        for (const auto& [key, uq] : self->state().computed_uq) {
            if (uq < 0.0) quenching_points++;
        }
        
        self->println("Final results: {} total points, {} quenching points ({:.1f}%)", 
                     total_computed, quenching_points, 
                     100.0 * quenching_points / total_computed);
        
        // Print refinement timing statistics
        self->println("Refinement timing statistics:");
        self->println("  Total refinement time: {:.6f}s ({} level operations)", 
                     self->state().total_refinement_time, self->state().refinement_time_measurements);
        
        if (self->state().refinement_time_measurements > 0) {
            double avg_refinement_time = self->state().total_refinement_time / self->state().refinement_time_measurements;
            self->println("  Average refinement time per level: {:.6f}s", avg_refinement_time);
        }
        
        // Print per-level timing breakdown
        for (size_t i = 0; i < self->state().level_refinement_times.size(); ++i) {
            self->println("  Level {} refinement time: {:.6f}s", i, self->state().level_refinement_times[i]);
        }
        
        // Trigger system termination check
        if (self->state().pending_jobs.empty() && self->state().active_jobs == 0) {
            self->println("Refinement complete and no active jobs - system should terminate");
        }
        return;
    }
    
    if (self->state().current_refinement_level > self->state().max_refinement_level) {
        self->println("=== REFINEMENT COMPLETE - Maximum level reached ===");
        
        // Clear boxes to trigger termination
        self->state().current_level_boxes.clear();
        
        if (self->state().pending_jobs.empty() && self->state().active_jobs == 0) {
            self->println("Max level reached and no active jobs - system should terminate");
        }
        return;
    }
    
    // Start processing new corner point jobs
    wake_idle_workers(self);
}

behavior server(stateful_actor<server_state>* self,
                uint16_t port,
                bool enable_bracket,
                bool enable_early_termination) {
  self->state().use_bracket = enable_bracket;
  self->state().use_early_termination = enable_early_termination;
  self->state().start_time  = std::chrono::high_resolution_clock::now();

  // Initialize system configuration with static parameters
  self->state().sys_config = system_config::create_static_config(enable_bracket, enable_early_termination);
  self->state().sys_config.finalize(); // Calculate derived parameters

  // Initialize computation configuration
  self->state().comp_config.enable_bracket_optimization = enable_bracket;
  self->state().comp_config.enable_early_termination = enable_early_termination;

  if (!self->system().middleman().publish(self, port)) {
    self->println("Failed to publish server on port {}", port);
    self->quit();
  }
  self->println("BASELINE Server up on port {}, bracket-enabled={}, early_termination={}", 
                port, enable_bracket, enable_early_termination);

  // Initialize spatial hash cell sizes based on coarse grid parameters
  self->state().spatial_mgr.initialize_cell_sizes(self->state().sys_config.theta_step, self->state().sys_config.xs_step);

  const size_t estimated_jobs = 2000;
  self->state().job_graph.reserve(estimated_jobs);
  self->state().scheduled_jobs.reserve(estimated_jobs);
  self->state().refinement_depth.reserve(estimated_jobs);
  
  // Reserve spatial hash grid capacity (estimate ~100-200 cells)
  self->state().spatial_mgr.grid.reserve(200);

  // Initialize box-based refinement system - TEST MODE: gamma=3 only, max 3 levels
  self->state().sys_config.gg_levels = {3, 2, 1, 0}; // Override to test only gamma=3
  self->println("=== BOX-BASED REFINEMENT TEST: gamma=3 only, {} workers ===", self->state().sys_config.actor_number);

  // Set max refinement level to 8 for testing
  self->state().max_refinement_level = 8;
  
  // Initialize 5x5 grid and create level 0 boxes
  initialize_coarse_grid_and_boxes(self);
  
  // Level 0 jobs are already created, just wake workers
  wake_idle_workers(self);
  
  self->println("=== BOX-BASED INITIALIZATION COMPLETE: {} pending jobs ===", 
                self->state().pending_jobs.size());

  // Spawn workers and send them initial jobs
  for (int i = 0; i < self->state().sys_config.actor_number; ++i) {
    auto worker = self->spawn(worker_actor, self, self->state().use_early_termination);
    self->state().active_workers.insert(worker);
    
    if (!self->state().pending_jobs.empty()) {
      job_t job = self->state().pending_jobs.front();
      self->state().pending_jobs.pop_front();

      if (should_apply_bracket_update(self, job)) {
        initialize_job_bracket(self, job);
      }

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
        
        // Bracket optimization happens at job creation time for new refinement jobs
        // For existing jobs being reassigned, no additional optimization needed
        
        self->anon_send(worker, std::move(job));
        ++self->state().active_jobs;
      } else {
        // No jobs available, keep worker idle
        self->state().idle_workers.insert(worker);
      }
    },

    [=](const std::string& msg, int gg, float theta, float xs, double uq, double usmin, double usmax, 
      double runtime_sec, actor worker) {
      try {
        if (msg != "result") return;

        // Store result in box-based system
        auto job_key = std::make_tuple(gg, theta, xs);
        self->state().computed_uq[job_key] = uq;
        
        --self->state().active_jobs;
        ++self->state().total_jobs_completed;
        
        // Update cache if uq < 0 (for bracket optimization)
        if (uq < 0.0) {
            self->state().uq_cache[{gg, theta}].emplace_back(xs, uq);
        }
        
        self->println("Job completed: gg={}, theta={:.3f}, xs={:.3f}, uq={:.6f}, bracket=[{:.3f},{:.3f}], runtime={:.3f}s", 
                      gg, theta, xs, uq, usmin, usmax, runtime_sec);

        // Assign next job to worker or make worker idle
        if (!self->state().pending_jobs.empty()) {
            job_t job = self->state().pending_jobs.front();
            self->state().pending_jobs.pop_front();
            ++self->state().active_jobs;
            
            // Apply bracket optimization for all eligible jobs
            if (should_apply_bracket_update(self, job)) {
                initialize_job_bracket(self, job);
            }
            
            self->anon_send(worker, std::move(job));
        } else {
            // No pending jobs - worker becomes idle
            self->state().idle_workers.insert(worker);
        }
        
        // Check if current refinement level is complete and advance if needed
        check_level_completion_and_advance(self);
        
        // Simple termination check: no pending jobs AND no active jobs = done
        if (!self->state().termination_initiated && 
            self->state().pending_jobs.empty() && 
            self->state().active_jobs == 0) {
            
            self->state().termination_initiated = true;
            
            self->println("=== TERMINATION: No pending jobs ({}) and no active jobs ({}) ===", 
                         self->state().pending_jobs.size(), self->state().active_jobs);
            
            // Send quit to all workers
            for (auto worker : self->state().active_workers) {
                self->anon_send(worker, "quit");
            }
            for (auto worker : self->state().idle_workers) {
                self->anon_send(worker, "quit");
            }
            
            int total_workers = self->state().active_workers.size() + self->state().idle_workers.size();
            self->println("Sent quit to {} workers. Waiting for confirmations...", total_workers);
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

void cleanup_completed_level_neighbors(stateful_actor<server_state>* self, int completed_level) {
  int cleaned_relationships = 0;
  
  // Collect all job keys from completed level to avoid iterator invalidation
  std::vector<std::tuple<int, float, float>> completed_level_jobs;
  for (auto& [job_key, job_node] : self->state().job_graph) {
    if (self->state().refinement_depth[job_key] == completed_level) {
      completed_level_jobs.push_back(job_key);
    }
  }
  
  // Clean up neighbor relationships for jobs from the completed level
  for (auto job_key : completed_level_jobs) {
    auto& job_node = self->state().job_graph[job_key];
    
    // Make a copy of neighbors to avoid iterator invalidation during erasure
    auto neighbors_copy = job_node.neighbors;
    
    for (auto neighbor_key : neighbors_copy) {
      int neighbor_depth = self->state().refinement_depth[neighbor_key];
      
      // Remove relationships between jobs at the same completed level
      // (they have all been processed and won't need to check each other again)
      if (neighbor_depth == completed_level) {
        job_node.neighbors.erase(neighbor_key);
        self->state().job_graph[neighbor_key].neighbors.erase(job_key);
        cleaned_relationships++;
      }
    }
  }
  
  if (cleaned_relationships > 0) {
    // self->println("NEIGHBOR CLEANUP: Removed {} obsolete neighbor relationships from level {}", 
    //              cleaned_relationships, completed_level);
  }
}

// MAIN
void caf_main(actor_system& system, const config& cfg) {
  scoped_actor self{system};

  self->println("Port: {} Host: {}", cfg.port, cfg.host);
  self->println("Server Mode: {}", cfg.server_mode);
  self->println("Early Termination (Sundials): {}", cfg.enable_early_termination);

  if (cfg.server_mode) {
    auto server_actor = system.spawn(server, cfg.port, cfg.enable_bracket, cfg.enable_early_termination);
    self->println("Server actor spawned");
  } else {
    auto remote_actor = system.spawn(remote, cfg.host, cfg.port, cfg.enable_early_termination);
    self->println("Remote actor spawned");
  }
}

CAF_MAIN(io::middleman, caf::id_block::my_project)