#include <cmath>
#include <chrono>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <optional>
#include <queue>
#include <ratio>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include <caf/io/all.hpp>

#include "config.h"
#include "system_config.h"  // Include the system configuration
#include "tuple_hash.h"  // Hash specializations for tuples
#include "job_structures.h"
#include "bracket_optimizer.h"
#include "dynamic_neighbor_provider.h"
#include "fitzhugh_nagumo.h"
#include "bisection_solver.h"
#include "wave_loader.h"

using namespace caf;
using namespace std::chrono_literals;

struct worker_state {
  // No need for complex initialization - sweep_main handles everything
};

struct remote_state {
  actor manager;
  bool use_early_termination;
};

behavior worker_actor(stateful_actor<worker_state>* self, actor manager, bool enable_early_termination);
behavior remote(stateful_actor<remote_state>* self, const std::string& hostname, uint16_t port, bool enable_early_termination);

// Use types from job_structures.h

// CAF serialization for job_t
template <class Inspector>
bool inspect(Inspector& f, job_t& x) {
  return f.object(x).fields(f.field("gg", x.gg),
                           f.field("theta", x.theta),
                           f.field("xs", x.xs),
                           f.field("n", x.n),
                           f.field("Ufile", x.Ufile),
                           f.field("Pfile", x.Pfile),
                           f.field("ufile", x.ufile),
                           f.field("pfile", x.pfile));
}

CAF_BEGIN_TYPE_ID_BLOCK(my_project, caf::first_custom_type_id)
CAF_ADD_TYPE_ID(my_project, (job_t))
CAF_END_TYPE_ID_BLOCK(my_project)

// --- SERVER STATE ---

struct server_state {

  // Use job_node from dynamic_neighbor_provider.h
  using job_node = bracket_optimizer::job_node;

  std::unordered_map<std::tuple<int, float, float>, job_node> job_graph;

  std::deque<job_t> pending_jobs;
  bracket_optimizer::uq_cache_t uq_cache;
  std::unordered_set<std::tuple<int, float, float>> scheduled_jobs;
  std::unordered_map<std::tuple<int, float, float>, int> refinement_depth;
  
  struct spatial_hash_manager {
    float theta_cell_size = 2.0f;  // Will be set based on coarse grid parameters
    float xs_cell_size = 10.0f;    // Will be set based on coarse grid parameters
    
    void initialize_cell_sizes(float theta_step, float xs_step) {
        theta_cell_size = theta_step / 5.0f;  
        xs_cell_size = xs_step / 3.5f;       
        
        // Ensure minimum reasonable cell sizes
        theta_cell_size = std::max(theta_cell_size, 0.5f);   // At least 0.5 degrees
        xs_cell_size = std::max(xs_cell_size, 2.0f);         // At least 2.0 xs units
    }
    
    // OPTIMIZATION: Use unordered_set instead of vector for O(1) removal
    std::unordered_map<std::tuple<int, int, int>, std::unordered_set<std::tuple<int, float, float>>> grid;
    
    std::tuple<int, int, int> get_cell_key(int gg, float theta, float xs) {
        int theta_cell = static_cast<int>(std::floor(theta / theta_cell_size));
        int xs_cell = static_cast<int>(std::floor(xs / xs_cell_size));
        return std::make_tuple(gg, theta_cell, xs_cell);
    }
    
    void add_job(const std::tuple<int, float, float>& job_key) {
        auto [gg, theta, xs] = job_key;
        auto cell_key = get_cell_key(gg, theta, xs);
        grid[cell_key].insert(job_key);  // O(1) insertion
    }
    
    void remove_job(const std::tuple<int, float, float>& job_key) {
        auto [gg, theta, xs] = job_key;
        auto cell_key = get_cell_key(gg, theta, xs);
        auto& cell_jobs = grid[cell_key];
        cell_jobs.erase(job_key);  // O(1) removal instead of O(n) erase-remove
        if (cell_jobs.empty()) {
            grid.erase(cell_key);
        }
    }
    
    std::vector<std::tuple<int, float, float>> find_neighbors(
        const std::tuple<int, float, float>& job_key, 
        float radius_theta, float radius_xs) {
        
        auto [gg, theta, xs] = job_key;
        std::vector<std::tuple<int, float, float>> neighbors;
        
        int cell_radius_theta = static_cast<int>(std::ceil(radius_theta / theta_cell_size)) + 1;
        int cell_radius_xs = static_cast<int>(std::ceil(radius_xs / xs_cell_size)) + 1;
        
        auto [center_gg, center_theta_cell, center_xs_cell] = get_cell_key(gg, theta, xs);
        
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
                        
                        // Use rectangular radius check (suitable for grid-based refinement)
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
  } spatial_mgr;

  bool use_bracket;
  bool use_early_termination;
  int  actor_number = 31;
  int  active_jobs = 0;

  // Computation configuration
  computation_config comp_config;
  
  // Bracket optimization components
  std::shared_ptr<bracket_optimizer::DynamicNeighborProvider> neighbor_provider;
  bracket_optimizer::BracketConfig bracket_config;

  // System configuration - single source of truth
  system_config conf;  // Pre-computed grid index lookups to avoid O(n) std::find_if operations
  std::unordered_map<float, int> theta_to_index;  // theta value -> grid index
  std::unordered_map<float, int> xs_to_index;     // xs value -> grid index
  std::vector<float> coarse_theta_values;         // Cached for reuse
  std::vector<float> coarse_xs_values;            // Cached for reuse

  // Worker tracking
  std::unordered_set<actor> active_workers;
  std::unordered_set<actor> idle_workers;
  
  // Job completion tracking
  int total_jobs_created = 0;
  int total_jobs_completed = 0;

  // Refinement timing tracking
  int total_refinement_jobs_created = 0;
  double total_refinement_time_ms = 0.0;  // Only job creation overhead
  int refinement_operations = 0;
  
  // Comprehensive refinement timing for dynamic approach
  std::chrono::high_resolution_clock::time_point first_refinement_start_time;
  bool refinement_timing_started = false;
  int initial_coarse_jobs = 0;
  int refinement_jobs_completed = 0;

  std::chrono::high_resolution_clock::time_point start_time;
};

// --- WORKER BEHAVIOR ---

behavior worker_actor(stateful_actor<worker_state>* self, actor manager, bool enable_early_termination) {
  return {
    [=](const exit_msg& msg) {
      if (msg.reason != exit_reason::kill) {
        // self->println("Worker {} received shutdown signal, exiting...", self->id());
      }
    },
    [=](job_t job) {
      auto start_time = std::chrono::high_resolution_clock::now();
      
      try {
        double usmin = -1000.0 / std::pow(2.0, job.n);
        double usmax = 0.0;
        
        // Use the sweep_main function which handles everything internally
        double result = sweep_main(
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
          enable_early_termination  // early termination flag
        );
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double runtime_sec = duration.count() / 1000.0;
        
        // Send result back to server
        self->mail("result", job.gg, job.theta, job.xs, result, 
                  usmin, usmax, runtime_sec, self).send(manager);
                  
      } catch (const std::exception& e) {
        self->println("Worker {} exception: {}", self->id(), e.what());
        
        // Send error result back
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double runtime_sec = duration.count() / 1000.0;
        
        self->mail("result", job.gg, job.theta, job.xs, 1.0, 
                  -1000.0, 0.0, runtime_sec, self).send(manager);
      }
    },

    [=](const std::string& msg) {
      if (msg == "quit") {
        // self->println("Worker {} shutting down...", self->id());
        self->mail("quit", self).send(manager);
        self->quit();
      }
    }
  };
}

// --- SERVER BEHAVIOR ---

void initialize_grid_lookups(stateful_actor<server_state>* self) {
  self->state().theta_to_index.clear();
  self->state().xs_to_index.clear();
  self->state().coarse_theta_values.clear();
  self->state().coarse_xs_values.clear();
  
  float theta_step = self->state().conf.theta_step;
  float xs_step = self->state().conf.xs_step;
  float max_theta = self->state().conf.max_theta;
  float xs_min = self->state().conf.xs_min;
  int NDNS = self->state().conf.NDNS;

  // Build theta lookups
  int theta_idx = 0;
  for (float theta = -max_theta; theta <= max_theta; theta += theta_step) {
    self->state().coarse_theta_values.push_back(theta);
    self->state().theta_to_index[theta] = theta_idx++;
  }
  
  // Build xs lookups
  for (int i = 0; i < NDNS; ++i) {
    float xs = xs_min + i * xs_step;
    self->state().coarse_xs_values.push_back(xs);
    self->state().xs_to_index[xs] = i;
  }
}

// Function to get grid index with O(1) lookup instead of O(n) std::find_if
int get_theta_index(stateful_actor<server_state>* self, float theta) {
  // Find closest theta value in coarse grid
  for (const auto& [grid_theta, index] : self->state().theta_to_index) {
    if (std::abs(grid_theta - theta) < 0.1f) {
      return index;
    }
  }
  return -1; 
}

int get_xs_index(stateful_actor<server_state>* self, float xs) {
  // Find closest xs value in coarse grid  
  for (const auto& [grid_xs, index] : self->state().xs_to_index) {
    if (std::abs(grid_xs - xs) < 1.0f) {
      return index;
    }
  }
  return -1; 
}

// Function to create neighbors for coarse grid jobs using static-like radius calculation
void create_coarse_grid_neighbors(stateful_actor<server_state>* self, 
                                 const std::tuple<int, float, float>& job_key) {
  
  auto [gg, theta, xs] = job_key;
  
  // Use static-like base radius calculation with system configuration parameters
  float base_theta_spacing = self->state().conf.theta_step;
  float base_xs_spacing = self->state().conf.xs_step;
  
  // Calculate base radius using configuration parameters
  float base_radius = std::min(base_theta_spacing, base_xs_spacing) * self->state().conf.base_radius_multiplier;
  
  // Apply coarse grid multiplier (depth 0 gets full base radius like static)
  float coarse_radius = base_radius * self->state().conf.coarse_grid_radius_multiplier;
  
  // Use euclidean radius for circular connectivity (like static)
  float euclidean_radius = coarse_radius * std::sqrt(2.0f);
  
  // Find neighbors using spatial hash
  auto neighbors = self->state().spatial_mgr.find_neighbors(job_key, 
                                                           euclidean_radius, 
                                                           euclidean_radius);
  
  // Connect to all found neighbors (naturally gives 8-connectivity for regular grid)
  for (auto& neighbor_key : neighbors) {
    auto [neighbor_gg, neighbor_theta, neighbor_xs] = neighbor_key;
    
    // Only connect to same gg and coarse grid jobs (depth = 0)
    if (neighbor_gg == gg && self->state().refinement_depth[neighbor_key] == 0) {
      // Calculate euclidean distance for circular connectivity
      float d_theta = neighbor_theta - theta;
      float d_xs = neighbor_xs - xs;
      float euclidean_distance = std::sqrt(d_theta * d_theta + d_xs * d_xs);
      
      // Connect if within euclidean radius (creates circular neighborhood like static)
      if (euclidean_distance <= euclidean_radius) {
        self->state().job_graph[job_key].neighbors.insert(neighbor_key);
        self->state().job_graph[neighbor_key].neighbors.insert(job_key);
      }
    }
  }
}

void schedule_job(stateful_actor<server_state>* self, job_t new_job, bool is_refinement = false) {
  auto key = std::make_tuple(new_job.gg, new_job.theta, new_job.xs);

  self->state().job_graph[key] = server_state::job_node{new_job, server_state::job_node::pending, 0.0, {}};
  
  self->state().total_jobs_created++;

  self->state().spatial_mgr.add_job(key);

  if (self->state().refinement_depth.find(key) == self->state().refinement_depth.end()) {
      self->state().refinement_depth[key] = 0; 
  }

  // Add to pending_jobs + scheduled_jobs (ready for workers + prevents duplicates)
  if (is_refinement) {
    self->state().pending_jobs.push_back(new_job);
  } else {
    self->state().pending_jobs.push_back(new_job);
  }
  self->state().scheduled_jobs.insert(key);
}


// Use the bracket optimizer from the separate module
void initialize_job_bracket(stateful_actor<server_state>* self, job_t& job) {
  if (self->state().use_bracket) {
    bracket_optimizer::initialize_bracket_from_cache(
      job,
      self->state().uq_cache,
      self->state().neighbor_provider,
      self->state().bracket_config
    );
  }
}

bool is_quenching_boundary(stateful_actor<server_state>* self,
                           const std::tuple<int, float, float>& job_key,
                           const std::tuple<int, float, float>& neighbor_key) {
  
  auto& job_node = self->state().job_graph[job_key];
  auto& neighbor_node = self->state().job_graph[neighbor_key];
  
  double uq1 = job_node.uq;
  double uq2 = neighbor_node.uq;
  
  bool job_quenches = (uq1 < 0.0);
  bool neighbor_quenches = (uq2 < 0.0);
  bool is_boundary = false;
  
  if (job_quenches != neighbor_quenches) {
    is_boundary = true;
  }
  
  // For both quenching: look for significant magnitude differences
  if (!is_boundary && job_quenches && neighbor_quenches) {
    int gg = std::get<0>(job_key); 
    double threshold = 1.0;  // Default threshold
    
    // Set gg-dependent thresholds based on observed uq ranges
    if (gg == 0) {
      threshold = 0.1;  
    } else if (gg == 1) {
      threshold = 0.2;  
    } else if (gg == 2) {
      threshold = 0.4;   
    } else if (gg == 3) {
      threshold = 1.0;  
    }
    
    if (std::abs(uq1 - uq2) > threshold) {
      is_boundary = true;
    }
  }

  // Remove this neighbor pair from both jobs' neighbor lists
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
  
  auto [gg1, theta1, xs1] = job_key;
  auto [gg2, theta2, xs2] = neighbor_key;

  // Only refine within the same gg value
  if (gg1 != gg2) {
    return;
  }

  // Calculate the distance between the two boundary points
  float d_theta = std::abs(theta1 - theta2);
  float d_xs = std::abs(xs1 - xs2);
  float boundary_distance = std::sqrt(d_theta * d_theta + d_xs * d_xs);
  
  // Calculate the radius of the new refinement point (half the distance)
  float refinement_radius = boundary_distance / 2.0f;
  
  // RADIUS-BASED STOPPING: Stop refining if the radius becomes smaller than h * factor
  // This ensures we don't refine beyond the finite difference discretization resolution
  float min_radius_threshold = self->state().conf.h * self->state().conf.min_refinement_radius_factor;
  if (refinement_radius < min_radius_threshold) {
    // self->println("REFINEMENT: Skipping - radius {:.6f} < threshold {:.6f} (h={:.6f})", 
    //               refinement_radius, min_radius_threshold, self->state().conf.h);
    return;
  }

  // Calculate midpoint between the two boundary jobs
  float center_theta = (theta1 + theta2) / 2.0f;
  float center_xs = (xs1 + xs2) / 2.0f;

  // Check bounds - use dynamic parameters
  if (center_theta < -self->state().conf.max_theta || center_theta > self->state().conf.max_theta ||
      center_xs < self->state().conf.xs_min || center_xs > self->state().conf.xs_max) {
    return;
  }
  
  auto new_job_key = std::make_tuple(gg1, center_theta, center_xs);
  if (self->state().job_graph.count(new_job_key) > 0) {
    return;
  }

  // Calculate the refinement depth for tracking purposes (still needed for logging/statistics)
  int depth1 = self->state().refinement_depth[job_key];
  int depth2 = self->state().refinement_depth[neighbor_key];
  int max_parent_depth = std::max(depth1, depth2);
  int new_depth = max_parent_depth + 1;
  
  // Calculate base spacing early for use in minimum distance calculations
  float base_theta_spacing = self->state().conf.theta_step;
  float base_xs_spacing = self->state().conf.xs_step;

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
  
  // Set refinement depth for new job BEFORE scheduling
  self->state().refinement_depth[new_job_key] = new_depth; 
  
  // Start timing the dynamic refinement process on first refinement job
  if (!self->state().refinement_timing_started) {
    self->state().first_refinement_start_time = std::chrono::high_resolution_clock::now();
    self->state().refinement_timing_started = true;
    self->state().initial_coarse_jobs = self->state().total_jobs_created; // Count before this refinement job
    // self->println("REFINEMENT TIMING: Starting dynamic refinement timing with {} initial coarse jobs", 
    //               self->state().initial_coarse_jobs);
  }
  
  // Schedule the job to add it to spatial hash and job graph (add to END of queue)
  schedule_job(self, new_job, true);  // true = is_refinement

  // DEPTH-BASED RADIUS REDUCTION: Apply static-like radius reduction based on refinement depth
  // This creates appropriate connectivity that shrinks with each refinement level
  
  // Calculate base radius using configuration parameters (same as coarse grid calculation)
  float base_radius = std::min(base_theta_spacing, base_xs_spacing) * self->state().conf.base_radius_multiplier;
  
  // Apply radius reduction factor based on refinement depth (60% retention per level like static)
  float radius_multiplier = std::pow(self->state().conf.radius_reduction_factor, new_depth);
  float search_radius = base_radius * radius_multiplier;
  
  // Use euclidean radius for circular connectivity (like static)
  float euclidean_radius = search_radius * std::sqrt(2.0f);
  
  // Debug output for radius calculation
  // self->println("Refinement depth {}: base_radius={:.1f}, search_radius={:.1f}, euclidean_radius={:.1f} (reduction={:.2f})", 
  //              new_depth, base_radius, search_radius, euclidean_radius, radius_multiplier);

  auto around = self->state().spatial_mgr.find_neighbors(new_job_key, euclidean_radius, euclidean_radius);
  
  // Connect to all jobs within the euclidean radius (creates circular neighborhood like static)
  int connections_made = 0;
  for (auto& existing_key : around) {
    auto [existing_gg, existing_theta, existing_xs] = existing_key;
    
    // Note: find_neighbors already filters by same gg, so no need to check again
    
    // Calculate euclidean distance in 2D parameter space (like static)
    float d_theta = existing_theta - center_theta;
    float d_xs = existing_xs - center_xs;
    float euclidean_distance = std::sqrt(d_theta * d_theta + d_xs * d_xs);
    
    // EUCLIDEAN RADIUS: Connect if within circular radius (creates proper circular neighborhood)
    if (euclidean_distance <= euclidean_radius) {
      self->state().job_graph[new_job_key].neighbors.insert(existing_key);
      self->state().job_graph[existing_key].neighbors.insert(new_job_key);
      connections_made++;
    }
  }
  
  // Record timing for this refinement operation
  auto refinement_end_time = std::chrono::high_resolution_clock::now();
  auto refinement_duration = std::chrono::duration_cast<std::chrono::microseconds>(refinement_end_time - refinement_start_time);
  double refinement_time_ms = refinement_duration.count() / 1000.0;
  
  self->state().total_refinement_time_ms += refinement_time_ms;
  self->state().total_refinement_jobs_created++;
  self->state().refinement_operations++;
  
  // self->println("REFINEMENT: Created job: gg={}, theta={}, xs={} (depth={}, euclidean_radius={:.3f}) with {} neighbors in {:.3f}ms", 
                // gg1, center_theta, center_xs, new_depth, euclidean_radius, connections_made, refinement_time_ms);
}

void create_all_coarse_grid_neighbors(stateful_actor<server_state>* self) {
  
  // Create all neighbor relationships for coarse grid jobs
  for (auto& [job_key, job_node] : self->state().job_graph) {    
    // Only process coarse grid jobs
    if (self->state().refinement_depth[job_key] == 0) {
      create_coarse_grid_neighbors(self, job_key);
    }
  }

  int total_connections = 0;
  for (auto& [job_key, job_node] : self->state().job_graph) {
    if (self->state().refinement_depth[job_key] == 0) {
      total_connections += job_node.neighbors.size();
    }
  }
  
}

behavior server(stateful_actor<server_state>* self,
                uint16_t port,
                bool enable_bracket,
                bool enable_early_termination) {
  // Initialize system configuration for dynamic mode
  self->state().conf = system_config::create_dynamic_config(enable_bracket, enable_early_termination);
  self->state().conf.finalize();  // Calculate derived parameters
  
  self->state().use_bracket = enable_bracket;
  self->state().use_early_termination = enable_early_termination;
  self->state().start_time  = std::chrono::high_resolution_clock::now();

  // Initialize computation configuration
  self->state().comp_config.enable_bracket_optimization = enable_bracket;
  self->state().comp_config.enable_early_termination = enable_early_termination;
  
  // Initialize bracket optimization components
  self->state().neighbor_provider = std::make_shared<bracket_optimizer::DynamicNeighborProvider>(&self->state().job_graph);
  self->state().bracket_config = bracket_optimizer::BracketConfig{};

  if (!self->system().middleman().publish(self, port)) {
    self->println("Failed to publish server on port {}", port);
    self->quit();
  }
      self->println("Server up on port {}, bracket-enabled={}, early-termination={}", 
                port, enable_bracket, enable_early_termination);

  // Initialize spatial hash cell sizes based on config parameters
  self->state().spatial_mgr.initialize_cell_sizes(self->state().conf.theta_step, self->state().conf.xs_step);

  // OPTIMIZATION: Initialize grid index lookups for O(1) access
  initialize_grid_lookups(self);

  std::vector<float> xs_values;
  for (int i = 0; i < self->state().conf.NDNS; ++i) {
    float xs = self->state().conf.xs_min + i * self->state().conf.xs_step;
    xs_values.push_back(xs);
  }

  std::vector<int> all_theta;
  for (float theta = -self->state().conf.max_theta; theta <= self->state().conf.max_theta; theta += self->state().conf.theta_step)
    all_theta.push_back(theta);

  // OPTIMIZATION: Pre-calculate job count and reserve hash map capacity to avoid rehashing
  // Number of jobs = gg_count * theta_count * xs_count
  std::vector<int> test_gg_values = {3, 2, 1, 0}; 
  size_t estimated_job_count = test_gg_values.size() * all_theta.size() * xs_values.size();
  size_t estimated_capacity = static_cast<size_t>(estimated_job_count * 1.5); // Extra capacity for refinement jobs
  
  self->state().job_graph.reserve(estimated_capacity);
  self->state().scheduled_jobs.reserve(estimated_capacity);
  self->state().refinement_depth.reserve(estimated_capacity);
  self->state().spatial_mgr.grid.reserve(estimated_capacity / 10);
  
  // Create jobs in reverse gg order (FIFO execution)
  for (int gg : test_gg_values) {
    std::string base = "/globalhome/tus210/HPC/quenchin_actor/";
    std::string waved = base + "waves/index_11";
    std::string crit = base + "waves/index_10";
    auto Ufile = waved + "/" + std::to_string(gg) + "/U";
    auto Pfile = waved + "/" + std::to_string(gg) + "/p";
    auto ufile = crit + "/" + std::to_string(gg) + "/U";
    auto pfile = crit + "/" + std::to_string(gg) + "/p";

    for (float xs : xs_values) {
      for (float theta : all_theta) {
        auto job_key = std::make_tuple(gg, theta, xs);
        self->state().scheduled_jobs.insert(job_key);

        job_t job{
          gg, theta, xs, 0,  // n=0 gives usmin=-1000, usmax=0
          Ufile, Pfile, ufile, pfile
        };
        schedule_job(self, job);
      }
    }
  }

  self->println("Enqueued {} jobs", self->state().pending_jobs.size());
  
  // Create complete coarse grid neighbor map after all jobs are scheduled
  create_all_coarse_grid_neighbors(self);

  for (int i = 0; i < self->state().actor_number; ++i) {
    auto worker = self->spawn(worker_actor, self, self->state().use_early_termination);
    if (!self->state().pending_jobs.empty()) {
      job_t job = self->state().pending_jobs.front();
      self->state().pending_jobs.pop_front();

      auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
      self->state().job_graph[job_key].status = server_state::job_node::running;
      self->anon_send(worker, std::move(job));
      ++self->state().active_jobs;
      
    } else {
      //we need to keep the worker alive for potential refinement jobs
      self->state().idle_workers.insert(worker);
      // self->println("Worker {} is now idle, waiting for potential jobs", worker.id());
      // self->anon_send(worker, "quit");
    }
  }


  return {
    [=](actor remote, int) {
      anon_mail(self->state().actor_number).send(remote);
    },

    [=](actor worker) {
      self->state().active_workers.insert(worker);
      
      if (!self->state().pending_jobs.empty()) {
        job_t job = self->state().pending_jobs.front();
        self->state().pending_jobs.pop_front();
        initialize_job_bracket(self, job);
        auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
        self->state().job_graph[job_key].status = server_state::job_node::running;

        ++self->state().active_jobs;
        self->anon_send(worker, std::move(job));
        
        // self->println("REGISTER: Sent job to worker {}, {} jobs remaining", 
        //               worker.id(), self->state().pending_jobs.size());
      } else {
        // No jobs available, but keep worker alive for potential refinement jobs
        self->state().idle_workers.insert(worker);
        // self->println("Worker {} is now idle, waiting for potential refinement jobs", worker.id());
      }
    },

    [=](const std::string& msg, int gg, float theta, float xs, double uq, 
      double usmin, double usmax, double runtime_sec, actor worker) {
      try {
        // self->println("SERVER DEBUG: Received message '{}' from worker {}", msg, worker.id());
        
        if (msg != "result") return;

        // self->println("SERVER DEBUG: Processing result for gg={}, theta={}, xs={}, uq={}, runtime={}s", 
        //               gg, theta, xs, uq, runtime_sec);
      
        auto job_key = std::make_tuple(gg, theta, xs);
      
        // self->println("RESULT RECEIVED: {} jobs pending before processing", 
        //               self->state().pending_jobs.size());
      
        --self->state().active_jobs;
        ++self->state().total_jobs_completed;
        
        // Track refinement job completion for dynamic timing
        auto depth_it = self->state().refinement_depth.find(job_key);
        int job_depth = (depth_it != self->state().refinement_depth.end()) ? depth_it->second : 0;
        if (job_depth > 0) { // This is a refinement job (not coarse grid)
            self->state().refinement_jobs_completed++;
        }

        // Update job_graph
        self->state().job_graph[job_key].status = server_state::job_node::done;
        self->state().job_graph[job_key].uq = uq;

        // Update cache if uq < 0
        if (uq < 0.0)
            self->state().uq_cache[{gg, theta}].emplace_back(xs, uq);

        int depth = 0;
        auto depth_lookup = self->state().refinement_depth.find(job_key);
        if (depth_lookup != self->state().refinement_depth.end())
          depth = depth_lookup->second;
        
        // Get neighbor count for this job
        int neighbor_count = self->state().job_graph[job_key].neighbors.size();

        self->println("Job completed: gg={}, theta={}, xs={}, result={}, usmin={}, usmax={}, runtime={}s, depth={}, neighbors={}", 
                      gg, theta, xs, uq, usmin, usmax, runtime_sec, depth, neighbor_count);

        // === IMMEDIATELY DISPATCH NEXT JOB TO KEEP WORKER BUSY ===
        if (!self->state().pending_jobs.empty()) { 
            job_t job = self->state().pending_jobs.front();
            self->state().pending_jobs.pop_front();
            ++self->state().active_jobs;
            
            self->println("Queue: {}, Active: {}, Completed: {}/{}, Remaining: {}", 
                          self->state().pending_jobs.size(), 
                          self->state().active_jobs,
                          self->state().total_jobs_completed,
                          self->state().total_jobs_created,
                          self->state().total_jobs_created - self->state().total_jobs_completed);

            initialize_job_bracket(self, job);

            auto next_job_key = std::make_tuple(job.gg, job.theta, job.xs);
            self->state().job_graph[next_job_key].status = server_state::job_node::running;
            self->anon_send(worker, std::move(job));
            
            // self->println("RESULT: Sent next job to worker {}, {} jobs now remaining", 
            //               worker.id(), self->state().pending_jobs.size());
        } else {
            // No pending jobs available right now - worker becomes idle
            self->state().idle_workers.insert(worker);
            // self->println("Worker {} now idle, {} active jobs remaining", worker.id(), self->state().active_jobs);
        }

        // Store number of pending jobs before boundary refinement (for tracking)
        size_t pending_before = self->state().pending_jobs.size();

        // === OPTIMIZED BOUNDARY PAIR PRE-FILTERING ===
        // Collect all boundary pairs first, then process them efficiently to avoid redundant checks
        std::unordered_set<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> boundary_pairs;
        std::unordered_set<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> processed_pairs;
        
        auto& node = self->state().job_graph[job_key];
        size_t total_neighbors = node.neighbors.size();
        
        // OPTIMIZATION: Reserve capacity based on neighbor count to avoid rehashing
        boundary_pairs.reserve(total_neighbors);
        processed_pairs.reserve(total_neighbors * 2); // Both directions
        
        // Diagnostic check for excessive neighbor count (should be rare with current approach)
        if (total_neighbors > 50) {
            self->println("WARNING: Job ({},{},{}) has {} neighbors - investigating...",
                          gg, theta, xs, total_neighbors);
            
            // Count neighbors by type for diagnosis
            int coarse_neighbors = 0, refinement_neighbors = 0;
            for (auto neighbor_key : node.neighbors) {
                int neighbor_depth = self->state().refinement_depth[neighbor_key];
                if (neighbor_depth == 0) coarse_neighbors++;
                else refinement_neighbors++;
            }
        }
        
        // Create a copy of neighbors to avoid iterator invalidation during is_quenching_boundary
        auto neighbors_copy = node.neighbors;
        
        for (auto neighbor_key : neighbors_copy) {
            auto& neighbor_node = self->state().job_graph[neighbor_key];

            // Only process neighbors that are already done
            if (neighbor_node.status == server_state::job_node::done) {
                // Check if this is a quenching boundary
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
        
        for (auto& boundary_pair : boundary_pairs) {
            auto job_key_pair = boundary_pair.first;
            auto neighbor_key_pair = boundary_pair.second;
            
            // Ensure we don't process the same pair twice
            if (processed_pairs.count(boundary_pair) == 0) {
                auto reverse_pair = std::make_pair(neighbor_key_pair, job_key_pair);
                processed_pairs.insert(boundary_pair);
                processed_pairs.insert(reverse_pair);  // Mark both directions as processed
                
                trigger_boundary_refinement(self, job_key_pair, neighbor_key_pair);
            }
        }

        // Check if new jobs were created and wake idle workers
        size_t new_jobs_created = self->state().pending_jobs.size() - pending_before;
        if (new_jobs_created > 0) {
            // self->println("REFINEMENT: Created {} new jobs, waking up idle workers", new_jobs_created);
            
            // Wake up idle workers for new refinement jobs
            auto idle_it = self->state().idle_workers.begin();
            while (idle_it != self->state().idle_workers.end() && !self->state().pending_jobs.empty()) {
                auto idle_worker = *idle_it;
                idle_it = self->state().idle_workers.erase(idle_it);
                
                job_t job = self->state().pending_jobs.front();
                self->state().pending_jobs.pop_front();
                ++self->state().active_jobs;
                
                initialize_job_bracket(self, job);
                
                auto wake_job_key = std::make_tuple(job.gg, job.theta, job.xs);
                self->state().job_graph[wake_job_key].status = server_state::job_node::running;
                self->anon_send(idle_worker, std::move(job));
                
                // self->println("WAKE: Sent refinement job to previously idle worker {}", idle_worker.id());
            }
        }

        if (self->state().pending_jobs.empty() && self->state().active_jobs == 0) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - self->state().start_time);
            self->println("=== ALL JOBS COMPLETED ===");
            self->println("Total runtime: {} seconds", duration.count());
            self->println("Jobs completed: {}/{}", self->state().total_jobs_completed, self->state().total_jobs_created);
            
            // Refinement timing summary
            self->println("=== REFINEMENT TIMING SUMMARY (DYNAMIC) ===");
            
            // Job creation overhead timing (original measurement)
            self->println("Refinement Job Creation Overhead:");
            self->println("  Total refinement operations: {}", self->state().refinement_operations);
            self->println("  Total refinement jobs created: {}", self->state().total_refinement_jobs_created);
            self->println("  Job creation time: {:.3f} ms", self->state().total_refinement_time_ms);
            if (self->state().refinement_operations > 0) {
                double avg_refinement_time = self->state().total_refinement_time_ms / self->state().refinement_operations;
                self->println("  Average job creation time: {:.3f} ms per operation", avg_refinement_time);
            }
            
            // Dynamic refinement process timing (comprehensive measurement)
            self->println("Dynamic Refinement Process Timing:");
            self->println("  Initial coarse jobs: {}", self->state().initial_coarse_jobs);
            self->println("  Refinement jobs completed: {}", self->state().refinement_jobs_completed);
            
            if (self->state().refinement_timing_started && self->state().refinement_jobs_completed > 0) {
                auto refinement_end_time = std::chrono::high_resolution_clock::now();
                auto total_refinement_duration = std::chrono::duration_cast<std::chrono::microseconds>(
                    refinement_end_time - self->state().first_refinement_start_time);
                double total_dynamic_refinement_time_ms = total_refinement_duration.count() / 1000.0;
                
                self->println("  Total dynamic refinement time: {:.2f} ms ({:.2f} seconds)", 
                              total_dynamic_refinement_time_ms, total_dynamic_refinement_time_ms / 1000.0);
                
                double avg_refinement_job_time = total_dynamic_refinement_time_ms / self->state().refinement_jobs_completed;
                self->println("  Average time per refinement job: {:.2f} ms", avg_refinement_job_time);
                
                // Dynamic refinement percentage
                double dynamic_refinement_percentage = (duration.count() > 0) ? 
                    (total_dynamic_refinement_time_ms / 1000.0) / duration.count() * 100.0 : 0.0;
                    
                self->println("Timing as Percentage of Total Runtime:");
                self->println("  Job creation overhead: {:.3f}%", 
                              (self->state().total_refinement_time_ms / 1000.0) / duration.count() * 100.0);
                self->println("  Complete dynamic refinement: {:.2f}%", dynamic_refinement_percentage);
            } else {
                self->println("  No refinement jobs were created or completed");
                self->println("Timing as Percentage of Total Runtime:");
                self->println("  Job creation overhead: {:.3f}%", 
                              (self->state().total_refinement_time_ms / 1000.0) / duration.count() * 100.0);
                self->println("  Complete dynamic refinement: 0.00%");
            }
            
            for (auto worker : self->state().active_workers) {
                anon_mail("quit").send(worker);
            }
            for (auto worker : self->state().idle_workers) {
                anon_mail("quit").send(worker);
            }
            
            self->quit();
        }
      
      } catch (const std::exception& e) {
        self->println("SERVER EXCEPTION in result handler: {}", e.what());
      } catch (...) {
        self->println("SERVER UNKNOWN EXCEPTION in result handler");
      }
    },

    [=](const std::string& msg, actor worker) {
      if (msg == "quit") {
        // self->println("SERVER DEBUG: Worker {} acknowledged quit", worker.id());
        self->state().active_workers.erase(worker);
        self->state().idle_workers.erase(worker);
        
        // Check if all workers have quit
        if (self->state().active_workers.empty() && self->state().idle_workers.empty()) {
          // self->println("SERVER: All workers have quit, server shutting down");
          self->quit();
        }
      }
    }
  };
}

// Remote behavior for client mode
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