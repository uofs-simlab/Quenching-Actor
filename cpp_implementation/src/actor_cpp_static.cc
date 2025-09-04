#include "config.h" 
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

constexpr float L = 2700.0f; 

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
      
    //   self->println("Starting Job gg={}, theta={}, xs={}, n={}, [u_min,u_max]=[{},{}], early_termination={}",
    //                 job.gg, job.theta, job.xs, job.n, usmin, usmax, self->state().use_early_termination);

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

struct server_state {

  struct job_node {
    job_t job;
    enum status_t { pending, running, done } status;
    double uq; 
    std::unordered_set<std::tuple<int, float, float>> neighbors; 
  };

  std::unordered_map<std::tuple<int, float, float>, job_node> job_graph;

  std::deque<job_t> pending_jobs;
  bracket_optimizer::uq_cache_t uq_cache;
  std::unordered_set<std::tuple<int, float, float>> scheduled_jobs;
  std::unordered_map<std::tuple<int, float, float>, int> refinement_depth;
  
  // Spatial hash manager for efficient neighbor queries
  spatial_hash_manager spatial_mgr;

  bool use_bracket;
  bool use_early_termination;
  int  actor_number = 31;
  int  active_jobs = 0;

  computation_config comp_config;
  
  bracket_optimizer::BracketConfig bracket_config;

  // Coarse grid parameters (for adaptive refinement)
  float theta_step = 20.0f;  
  float xs_step = 50.0f; 
  float max_theta = 80.0f;
  float xs_min = 0.05f;
  float xs_max = 1000.0f;
  int NDNS = 10;

  std::unordered_map<float, int> theta_to_index;
  std::unordered_map<float, int> xs_to_index;

  // Worker tracking
  std::unordered_set<actor> active_workers;
  std::unordered_set<actor> idle_workers;
  
  // Job completion tracking
  int total_jobs_created = 0;
  int total_jobs_completed = 0;

  
  std::vector<int> gg_levels = {3, 2, 1, 0};  
  int global_refinement_level = 0;                // Current refinement level for ALL gg values (0=coarse, 1=first refinement, etc.)
  std::unordered_map<int, bool> gg_level_processing;        // Whether a gg level is currently processing jobs at current refinement level
  std::unordered_map<int, bool> gg_completely_done;         // Whether gg level is completely finished (no more refinement possible)
  
  // Job tracking per gg and refinement level
  std::unordered_map<std::pair<int, int>, int> gg_level_jobs_created;    
  std::unordered_map<std::pair<int, int>, int> gg_level_jobs_completed;  

  const int max_refinement_level = 8;  
  std::chrono::high_resolution_clock::time_point start_time;
};


void initialize_job_bracket(stateful_actor<server_state>* self, job_t& job) {
  if (self->state().use_bracket) {
    bracket_optimizer::initialize_bracket_from_cache(
      job,
      self->state().uq_cache,
      nullptr,  // No neighbor provider for static version
      self->state().bracket_config
    );
  }
}

// Helper function to determine if bracket updates should be applied
bool should_apply_bracket_update(stateful_actor<server_state>* self, const job_t& job) {
  if (!self->state().use_bracket) return false;
  
  // Get the first gg level in the sequence
  int first_gg = self->state().gg_levels[0];
  
  
  auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
  auto depth_it = self->state().refinement_depth.find(job_key);
  int depth = (depth_it != self->state().refinement_depth.end()) ? depth_it->second : 0;
  
  // Skip bracket updates only for coarse grid jobs (depth=0) of first gg
  if (job.gg == first_gg && depth == 0) {
    return false; 
  }
  return true; 
}

bool is_quenching_boundary(stateful_actor<server_state>* self, const std::tuple<int, float, float>& job_key, const std::tuple<int, float, float>& neighbor_key);
void trigger_boundary_refinement(stateful_actor<server_state>* self, const std::tuple<int, float, float>& job_key, const std::tuple<int, float, float>& neighbor_key);
void schedule_gg_level_jobs(stateful_actor<server_state>* self, int gg);
void check_gg_completion_and_refine(stateful_actor<server_state>* self);
void create_all_coarse_grid_neighbors(stateful_actor<server_state>* self);
void wake_idle_workers(stateful_actor<server_state>* self);


void initialize_grid_lookups(stateful_actor<server_state>* self) {
  self->state().theta_to_index.clear();
  int theta_idx = 0;
  for (float theta = -self->state().max_theta; theta <= self->state().max_theta; theta += self->state().theta_step) {
    self->state().theta_to_index[theta] = theta_idx++;
  }
  
  self->state().xs_to_index.clear();
  for (int i = 0; i < self->state().NDNS; ++i) {
    float xs = self->state().xs_min + i * self->state().xs_step;
    self->state().xs_to_index[xs] = i;
  }
}

int get_theta_index(stateful_actor<server_state>* self, float theta) {
  auto it = self->state().theta_to_index.find(theta);
  return (it != self->state().theta_to_index.end()) ? it->second : -1;
}

int get_xs_index(stateful_actor<server_state>* self, float xs) {
  auto it = self->state().xs_to_index.find(xs);
  return (it != self->state().xs_to_index.end()) ? it->second : -1;
}

void create_coarse_grid_neighbors(stateful_actor<server_state>* self, 
                                 const std::tuple<int, float, float>& job_key) {
  
  auto [gg, thetaC, xsC] = job_key;
  
  for (auto& [keyE, nodeE] : self->state().job_graph) {
      if (keyE == job_key) continue; // skip self
      if (std::get<0>(keyE) != gg) continue; // only same gg
      
      // Only connect to other coarse jobs (refinement_depth == 0)
      if (self->state().refinement_depth[keyE] != 0) continue;

      auto [_, thetaE, xsE] = keyE;

      int theta_idx = get_theta_index(self, thetaC);
      int theta_idx_E = get_theta_index(self, thetaE);
      int xs_idx = get_xs_index(self, xsC);
      int xs_idx_E = get_xs_index(self, xsE);
      
      // Check if both jobs are on the coarse grid
      if (theta_idx != -1 && theta_idx_E != -1 && xs_idx != -1 && xs_idx_E != -1) {
          
          // EXPLICIT 8-NEIGHBOR CONNECTIVITY (includes all diagonals)
          // Calculate grid distance for both dimensions
          int delta_theta = std::abs(theta_idx - theta_idx_E);
          int delta_xs = std::abs(xs_idx - xs_idx_E);
          
          // Connect to all 8 surrounding neighbors:
          bool is_horizontal = (delta_theta == 1 && delta_xs == 0);
          bool is_vertical = (delta_theta == 0 && delta_xs == 1);
          bool is_diagonal = (delta_theta == 1 && delta_xs == 1);
          
          if (is_horizontal || is_vertical || is_diagonal) {
              // Add bidirectional connection
              self->state().job_graph[job_key].neighbors.insert(keyE);
              self->state().job_graph[keyE].neighbors.insert(job_key);
          }
      }
  }
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
//   self->println("Scheduled job: gg={}, theta={}, xs={:.3f}, depth={}, type={}, total_created={}", 
//                 new_job.gg, new_job.theta, new_job.xs, depth, 
//                 is_refinement ? "refinement" : "coarse", self->state().total_jobs_created);
}
void check_gg_completion_and_refine(stateful_actor<server_state>* self) {
  
  bool all_gg_completed_current_level = true;
  int current_level = self->state().global_refinement_level;
  
  for (int gg : self->state().gg_levels) {
    // Skip if completely done
    if (self->state().gg_completely_done[gg]) continue;
    
    auto level_key = std::make_pair(gg, current_level);
    int level_created = self->state().gg_level_jobs_created[level_key];
    int level_completed = self->state().gg_level_jobs_completed[level_key];
    
    // DEBUG: Print status for each gg level
    // self->println("=== GG {} LEVEL {} STATUS: jobs={}/{}, processing={}, done={} ===", 
    //               gg, current_level, level_completed, level_created,
    //               self->state().gg_level_processing[gg], self->state().gg_completely_done[gg]);
    
    if (level_completed < level_created || level_created == 0) {
      all_gg_completed_current_level = false;
    }
  }
  
  // If ALL gg levels completed current level, advance ALL to next level
  if (all_gg_completed_current_level) {
    // self->println("=== ALL GG LEVELS COMPLETED LEVEL {} - ADVANCING ALL TO LEVEL {} ===", 
    //               current_level, current_level + 1);
    
    // Check if we can advance to next refinement level
    int next_level = current_level + 1;
    if (next_level > self->state().max_refinement_level) {
      // All levels reached maximum refinement
    //   self->println("=== ALL GG LEVELS REACHED MAXIMUM REFINEMENT LEVEL {} ===", self->state().max_refinement_level);
      
      // Mark all as completely done
      for (int gg : self->state().gg_levels) {
        self->state().gg_completely_done[gg] = true;
      }
    } else {
      // Advance global refinement level
      self->state().global_refinement_level = next_level;
      
      // Create refinement jobs for ALL gg levels at the new level
      int total_refinement_jobs = 0;
      
      for (int gg : self->state().gg_levels) {
        if (!self->state().gg_completely_done[gg]) {
          // Mark as processing next level
          self->state().gg_level_processing[gg] = true;
          
          // Find boundary pairs for this gg at the completed level
          std::set<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> boundary_pairs;
          std::set<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> processed_pairs;
          
          // Scan all completed jobs in this gg level for boundaries
          for (auto& [job_key, job_node] : self->state().job_graph) {
            auto [job_gg, job_theta, job_xs] = job_key;
            if (job_gg != gg || job_node.status != server_state::job_node::done) continue;
            
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
                
                if (max_parent_depth >= self->state().max_refinement_level) {
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
          
        //   self->println("=== GG {} advanced to level {} with {} jobs created ===", 
        //                 gg, next_level, gg_refinement_jobs);
          
          // If no jobs created for this gg, mark as completely done
          if (gg_refinement_jobs == 0) {
            self->state().gg_completely_done[gg] = true;
            // self->println("=== GG {} COMPLETELY FINISHED - no boundary refinement possible ===", gg);
          }
          
          self->state().gg_level_processing[gg] = false;
        }
      }
      
      if (total_refinement_jobs > 0) {
        // self->println("=== SYNCHRONIZED REFINEMENT: Advanced ALL gg levels to level {} with {} total jobs ===", 
                    //   next_level, total_refinement_jobs);
        
        // Wake workers for new jobs
        wake_idle_workers(self);
      } else {
        // No refinement jobs created for any gg - all are done
        // self->println("=== ALL GG LEVELS COMPLETELY FINISHED - no more refinement possible ===");
        for (int gg : self->state().gg_levels) {
          self->state().gg_completely_done[gg] = true;
        }
      }
    }
  }
  
  // Check if ALL gg levels are done
  bool all_done = true;
  for (int gg : self->state().gg_levels) {
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
    
    self->quit();
  }
  wake_idle_workers(self);
}

void wake_idle_workers(stateful_actor<server_state>* self) {
  auto idle_it = self->state().idle_workers.begin();
  while (idle_it != self->state().idle_workers.end() && !self->state().pending_jobs.empty()) {
    auto idle_worker = *idle_it;
    idle_it = self->state().idle_workers.erase(idle_it);
    
    job_t job = self->state().pending_jobs.front();
    self->state().pending_jobs.pop_front();
    ++self->state().active_jobs;
    
    if (should_apply_bracket_update(self, job))
      initialize_job_bracket(self, job);
    
    auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
    self->state().job_graph[job_key].status = server_state::job_node::running;
    self->anon_send(idle_worker, std::move(job));
  }
}

// Schedule jobs for a specific gg level
void schedule_gg_level_jobs(stateful_actor<server_state>* self, int gg) {
  float max_theta = self->state().max_theta;
  int theta_step = self->state().theta_step;
  const float xs_min = self->state().xs_min;
  const float xs_max = self->state().xs_max;
  int NDNS = self->state().NDNS;
  
  std::vector<float> xs_values;
  for (int i = 0; i < NDNS; ++i) {
    float xs = xs_min + i * self->state().xs_step;  // Use pre-calculated xs_step
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
  
//   self->println("Scheduled {} coarse grid jobs for gg level {}", jobs_for_gg, gg);
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
  
  // For both quenching: look for significant magnitude differences
  if (!is_boundary && job_quenches && neighbor_quenches) {
    // Dynamic threshold based on gg value - different gg values have different uq ranges
    int gg = std::get<0>(job_key);  // Both jobs should have same gg value
    double threshold = 1.0; 
    
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
    //   self->println("*** QUENCH MAGNITUDE BOUNDARY: ({},{},{}) uq={} vs ({},{},{}) uq={} diff={} (threshold={}) ***", 
    //                 std::get<0>(job_key), std::get<1>(job_key), std::get<2>(job_key), uq1,
    //                 std::get<0>(neighbor_key), std::get<1>(neighbor_key), std::get<2>(neighbor_key), uq2,
    //                 std::abs(uq1 - uq2), threshold);
      is_boundary = true;
    }
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
  int max_depth = self->state().max_refinement_level;  // Use dynamic max refinement level
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
  if (center_theta < -self->state().max_theta || center_theta > self->state().max_theta ||
      center_xs < self->state().xs_min || center_xs > self->state().xs_max) {
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
  
  // Schedule the job to add it to spatial hash and job graph (add to END of queue)
  schedule_job(self, new_job, true);  // true = is_refinement

  // Track refinement jobs per gg level (use global refinement level)
  int current_level = self->state().global_refinement_level;
  auto level_key = std::make_pair(gg1, current_level);
  self->state().gg_level_jobs_created[level_key]++;

  // OPTIMIZATION: ADAPTIVE NEIGHBOR RADIUS for better distribution
  // Use depth-adaptive radius that gradually decreases but maintains connectivity
  // This balances refinement precision with neighbor connectivity
  float base_theta_spacing = self->state().theta_step;  
  float base_xs_spacing = self->state().xs_step;        
  
  // Adaptive radius that decreases with refinement depth
  // depth_factor starts at 1.0 (depth 0) and gradually decreases
  float depth_factor = 1.0f / (1.0f + 0.1f * new_depth);  // Gradual decrease
  
  // Apply adaptive scaling with minimum connectivity guarantees
  float adaptive_radius_factor = 0.6f * depth_factor;  // Base 60% scaled by depth
  float radius_theta = base_theta_spacing * adaptive_radius_factor;
  float radius_xs = base_xs_spacing * adaptive_radius_factor;
  
  // Add a small tolerance factor to ensure we catch boundary neighbors
  // This is especially important for floating-point precision issues
  float tolerance_factor = 1.2f;  // 20% extra radius for safety
  radius_theta *= tolerance_factor;
  radius_xs *= tolerance_factor;
  
  
  // Find all existing jobs within the adaptive radius
  auto around = self->state().spatial_mgr.find_neighbors(new_job_key, radius_theta, radius_xs);
  
  // Connect to all jobs within the fixed radius
  // This naturally includes horizontal, vertical, and diagonal connections
  int connections_made = 0;
  for (auto& existing_key : around) {
    auto [existing_gg, existing_theta, existing_xs] = existing_key;
    
    // Note: find_neighbors already filters by same gg, so no need to check again
    
    // Calculate actual distances for both dimensions
    float d_theta = std::abs(existing_theta - center_theta);
    float d_xs = std::abs(existing_xs - center_xs);
    
    // ADAPTIVE RECTANGULAR RADIUS: Connect if within BOTH radii
    // This creates natural 8-directional connectivity (horizontal, vertical, diagonal)
    if (d_theta <= radius_theta && d_xs <= radius_xs) {
      self->state().job_graph[new_job_key].neighbors.insert(existing_key);
      self->state().job_graph[existing_key].neighbors.insert(new_job_key);
      connections_made++;
    }
  }
  
//   self->println("REFINEMENT: Created job (gg={} global_level {} job #{}): gg={}, theta={}, xs={} (depth={}) with {} neighbors", 
//                 gg1, current_level, self->state().gg_level_jobs_created[level_key], gg1, center_theta, center_xs, new_depth, connections_made);
}

// Function to create complete coarse grid neighbor map after all jobs are registered
void create_all_coarse_grid_neighbors(stateful_actor<server_state>* self) {
  
//   self->println("Creating coarse grid neighbor map for {} total jobs...", self->state().job_graph.size());
  
  int neighbors_created = 0;
  
  // Create all neighbor relationships for coarse grid jobs
  for (auto& [job_key, job_node] : self->state().job_graph) {    
    if (self->state().refinement_depth[job_key] == 0) {
      int initial_neighbors = job_node.neighbors.size();
      create_coarse_grid_neighbors(self, job_key);
      int final_neighbors = job_node.neighbors.size();
      neighbors_created += (final_neighbors - initial_neighbors);
      
      // Debug: Print neighbor count for each coarse job
      // auto [gg, theta, xs] = job_key;
    //   self->println("Coarse job gg={}, theta={}, xs={:.1f} has {} neighbors", 
    //                 gg, theta, xs, final_neighbors);
    }
  }
  
//   self->println("Coarse grid neighbor map complete. Total new neighbor connections: {}", neighbors_created);
}

behavior server(stateful_actor<server_state>* self,
                uint16_t port,
                bool enable_bracket,
                bool enable_early_termination) {
  self->state().use_bracket = enable_bracket;
  self->state().use_early_termination = enable_early_termination;
  self->state().start_time  = std::chrono::high_resolution_clock::now();

  // Initialize computation configuration
  self->state().comp_config.enable_bracket_optimization = enable_bracket;
  self->state().comp_config.enable_early_termination = enable_early_termination;
  
  // Initialize bracket optimization components
  self->state().bracket_config = bracket_optimizer::BracketConfig{};

  if (!self->system().middleman().publish(self, port)) {
    self->println("Failed to publish server on port {}", port);
    self->quit();
  }
  self->println("BASELINE Server up on port {}, bracket-enabled={}, early_termination={}", 
                port, enable_bracket, enable_early_termination);

  self->state().xs_step = (self->state().xs_max - self->state().xs_min) / (self->state().NDNS - 1);

  // Initialize spatial hash cell sizes based on coarse grid parameters
  self->state().spatial_mgr.initialize_cell_sizes(self->state().theta_step, self->state().xs_step);

  // OPTIMIZATION: Initialize grid lookup maps for O(1) index access
  initialize_grid_lookups(self);

  const size_t estimated_jobs = 2000;
  self->state().job_graph.reserve(estimated_jobs);
  self->state().scheduled_jobs.reserve(estimated_jobs);
  self->state().refinement_depth.reserve(estimated_jobs);
  
  // Reserve spatial hash grid capacity (estimate ~100-200 cells)
  self->state().spatial_mgr.grid.reserve(200);

  self->println("=== SYNCHRONIZED GG REFINEMENT: Initializing {} gg levels at global level {} ===", 
                self->state().gg_levels.size(), self->state().global_refinement_level);
  
  for (int gg : self->state().gg_levels) {
    self->state().gg_level_processing[gg] = false;
    self->state().gg_completely_done[gg] = false;
    
    // Schedule coarse grid jobs for this gg
    self->println("Scheduling coarse grid jobs for gg level {}...", gg);
    schedule_gg_level_jobs(self, gg);
    
    // Job counters already initialized in schedule_gg_level_jobs
    auto level_0_key = std::make_pair(gg, 0);
    auto gg_coarse_jobs = self->state().gg_level_jobs_created[level_0_key];
    
    self->println("GG {} initialized: global level 0 with {} coarse jobs", gg, gg_coarse_jobs);
  }

  // Create complete coarse grid neighbor map for ALL scheduled jobs
  create_all_coarse_grid_neighbors(self);

  self->println("=== INITIALIZATION COMPLETE: {} total jobs from {} gg levels ===", 
                self->state().pending_jobs.size(), self->state().gg_levels.size());

  // Spawn workers and send them initial jobs
  for (int i = 0; i < self->state().actor_number; ++i) {
    auto worker = self->spawn(worker_actor, self, self->state().use_early_termination);
    self->state().active_workers.insert(worker);
    
    if (!self->state().pending_jobs.empty()) {
      job_t job = self->state().pending_jobs.front();
      self->state().pending_jobs.pop_front();

      if (should_apply_bracket_update(self, job))
        initialize_job_bracket(self, job);

      auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
      self->state().job_graph[job_key].status = server_state::job_node::running;
      self->anon_send(worker, std::move(job));
      ++self->state().active_jobs;
    } else {
      self->state().idle_workers.insert(worker);
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
        
        if (should_apply_bracket_update(self, job))
          initialize_job_bracket(self, job);
          
        auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
        self->state().job_graph[job_key].status = server_state::job_node::running;
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

        auto job_key = std::make_tuple(gg, theta, xs);
        --self->state().active_jobs;
        ++self->state().total_jobs_completed;
        auto depth_it = self->state().refinement_depth.find(job_key);
        int depth = (depth_it != self->state().refinement_depth.end()) ? depth_it->second : 0;
        
        int tracking_level = self->state().global_refinement_level;
        
        auto level_key = std::make_pair(gg, tracking_level);
        ++self->state().gg_level_jobs_completed[level_key];
        
        // self->println("*** JOB COMPLETED for gg={} global_level={}: {}/{} jobs done ***", 
        //               gg, tracking_level, 
        //               self->state().gg_level_jobs_completed[level_key], 
        //               self->state().gg_level_jobs_created[level_key]);

        // Update job_graph
        self->state().job_graph[job_key].status = server_state::job_node::done;
        self->state().job_graph[job_key].uq = uq;

        // Update cache if uq < 0
        if (uq < 0.0)
            self->state().uq_cache[{gg, theta}].emplace_back(xs, uq);

        int neighbor_count = self->state().job_graph[job_key].neighbors.size();
        
        self->println("Job completed: gg={}, theta={}, xs={}, result={}, usmin={}, usmax={}, runtime={}s, depth={}, neighbors={}", 
                      gg, theta, xs, uq, usmin, usmax, runtime_sec, depth, neighbor_count);

        // === AUTOMATIC RESULT DEDUCTION ===
        // COMMENTED OUT: Auto-deduction creates inconsistency with boundary detection
        // Jobs that are deduced (not actually run) don't trigger boundary refinement
        // This causes different refinement patterns - disabling to ensure fair comparison
        /*
        if (uq >= 0.0) {  // Job did not quench (result = 1.0)
            for (auto it = self->state().pending_jobs.begin(); it != self->state().pending_jobs.end(); ) {
                const job_t& pending_job = *it;
                
                // Same gg and theta, but smaller xs
                if (pending_job.gg == gg && 
                    std::abs(pending_job.theta - theta) < 1e-4f && 
                    pending_job.xs < xs) {
                    
                    auto pending_key = std::make_tuple(pending_job.gg, pending_job.theta, pending_job.xs);
                    
                    // Mark as completed in job_graph
                    self->state().job_graph[pending_key].status = server_state::job_node::done;
                    self->state().job_graph[pending_key].uq = 1.0;  // Non-quenched result
                    
                    // Get depth and neighbor count for deduced job
                    auto pending_depth_it = self->state().refinement_depth.find(pending_key);
                    int pending_depth = (pending_depth_it != self->state().refinement_depth.end()) ? pending_depth_it->second : 0;
                    int pending_neighbor_count = self->state().job_graph[pending_key].neighbors.size();
                    
                    // Print the result like a normal job completion
                    self->println("Job completed: gg={}, theta={}, xs={}, result={}, usmin={}, usmax={}, runtime={}s, depth={}, neighbors={}", 
                                  pending_job.gg, pending_job.theta, pending_job.xs, 1.0, usmin, usmax, 0.0, pending_depth, pending_neighbor_count);

                    // Remove from scheduled jobs
                    self->state().scheduled_jobs.erase(pending_key);
                    
                    // Increment completion counters
                    ++self->state().total_jobs_completed;
                    
                    // Update the (gg, refinement_level) specific counters
                    auto pending_level_key = std::make_pair(pending_job.gg, pending_depth);
                    ++self->state().gg_level_jobs_completed[pending_level_key];
                    
                    // Remove from pending queue
                    it = self->state().pending_jobs.erase(it);
                } else {
                    ++it;
                }
            }
        }
        */

        if (!self->state().pending_jobs.empty()) { 
            job_t job = self->state().pending_jobs.front();
            self->state().pending_jobs.pop_front();
            ++self->state().active_jobs;
            
            // self->println("Queue: {}, Active: {}, Completed: {}/{}, Remaining: {}", 
            //               self->state().pending_jobs.size(), 
            //               self->state().active_jobs,
            //               self->state().total_jobs_completed,
            //               self->state().total_jobs_created,
            //               self->state().total_jobs_created - self->state().total_jobs_completed);

            if (should_apply_bracket_update(self, job))
                initialize_job_bracket(self, job);
                
            auto next_job_key = std::make_tuple(job.gg, job.theta, job.xs);
            self->state().job_graph[next_job_key].status = server_state::job_node::running;
            self->anon_send(worker, std::move(job));
        } else {
            // No pending jobs available - worker becomes idle
            self->state().idle_workers.insert(worker);
            
            // DEBUG: Print when worker becomes idle and check if all workers are idle
            int total_workers = self->state().active_workers.size() + self->state().idle_workers.size();
            // self->println("Worker became idle. Idle workers: {}/{}, Pending jobs: {}", 
            //               self->state().idle_workers.size(), total_workers, self->state().pending_jobs.size());
        }

        // Check if any gg level can advance to next refinement level
        check_gg_completion_and_refine(self);
        
        if (self->state().idle_workers.size() == self->state().active_workers.size() + self->state().idle_workers.size() && 
            self->state().pending_jobs.empty()) {
                bool any_incomplete = false;
                for (int gg : self->state().gg_levels) {
                    if (!self->state().gg_completely_done[gg]) {
                        any_incomplete = true;
                        break;
                    }
                }
                if (any_incomplete) {
                    self->println("*** WARNING: All workers idle but some gg levels incomplete! ***");
                    // Force a completion check to see if we can generate more jobs
                    check_gg_completion_and_refine(self);
                }
        }

        // Check completion and refinement for all gg levels
        check_gg_completion_and_refine(self);
      
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
        
        // Check if all workers have quit
        if (self->state().active_workers.empty() && self->state().idle_workers.empty()) {
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
    auto server_actor = system.spawn(server, cfg.port, cfg.enable_bracket, cfg.enable_early_termination);
    self->println("Server actor spawned");
  } else {
    auto remote_actor = system.spawn(remote, cfg.host, cfg.port, cfg.enable_early_termination);
    self->println("Remote actor spawned");
  }
}

CAF_MAIN(io::middleman, caf::id_block::my_project)