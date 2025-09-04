#include "config.h" 
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


constexpr float L = 2700.0f;  // System length parameter used throughout the simulation

// Hash functions for tuple types used as keys
namespace std {
  template<>
  struct hash<std::tuple<int, float, float>> {
    std::size_t operator()(const std::tuple<int, float, float>& t) const {
      auto h1 = std::hash<int>{}(std::get<0>(t));
      auto h2 = std::hash<float>{}(std::get<1>(t));
      auto h3 = std::hash<float>{}(std::get<2>(t));
      // FNV-1a inspired hash combining with golden ratio for better distribution
      std::size_t result = h1;
      result ^= h2 + 0x9e3779b9 + (result << 6) + (result >> 2);
      result ^= h3 + 0x9e3779b9 + (result << 6) + (result >> 2);
      return result;
    }
  };
  
  template<>
  struct hash<std::tuple<int, float>> {
    std::size_t operator()(const std::tuple<int, float>& t) const {
      auto h1 = std::hash<int>{}(std::get<0>(t));
      auto h2 = std::hash<float>{}(std::get<1>(t));
      // FNV-1a inspired hash combining
      std::size_t result = h1;
      result ^= h2 + 0x9e3779b9 + (result << 6) + (result >> 2);
      return result;
    }
  };
  
  template<>
  struct hash<std::tuple<int, int, int>> {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
      auto h1 = std::hash<int>{}(std::get<0>(t));
      auto h2 = std::hash<int>{}(std::get<1>(t));
      auto h3 = std::hash<int>{}(std::get<2>(t));
      // FNV-1a inspired hash combining
      std::size_t result = h1;
      result ^= h2 + 0x9e3779b9 + (result << 6) + (result >> 2);
      result ^= h3 + 0x9e3779b9 + (result << 6) + (result >> 2);
      return result;
    }
  };
  
  template<>
  struct hash<std::pair<int, int>> {
    std::size_t operator()(const std::pair<int, int>& p) const {
      auto h1 = std::hash<int>{}(p.first);
      auto h2 = std::hash<int>{}(p.second);
      // FNV-1a inspired hash combining
      std::size_t result = h1;
      result ^= h2 + 0x9e3779b9 + (result << 6) + (result >> 2);
      return result;
    }
  };
} 

//JOB & CACHE TYPES
struct job_t {
  int      gg;
  float    theta;
  float    xs;
  int      n;           // Power for bracket: usmin = -1000/(2^n), usmax = 0
  std::string Ufile;
  std::string Pfile;
  std::string ufile;
  std::string pfile;

  template <class Inspector>
  friend bool inspect(Inspector& f, job_t& x) {
    return f.object(x).fields(
      f.field("gg",     x.gg),
      f.field("theta",  x.theta),
      f.field("xs",     x.xs),
      f.field("n",      x.n),
      f.field("Ufile",  x.Ufile),
      f.field("Pfile",  x.Pfile),
      f.field("ufile",  x.ufile),
      f.field("pfile",  x.pfile)
      // f.field("savedir",x.savedir)
    );
  }
};

struct result_entry {
  double xs;
  double uq;
};

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
      
      // Compute usmin and usmax from n for display and Julia call
      double usmin = -1000.0 / std::pow(2.0, job.n);
      double usmax = 0.0;
      
      self->println("Starting Job gg={}, theta={}, xs={}, n={}, [u_min,u_max]=[{},{}], early_termination={}",
                    job.gg, job.theta, job.xs, job.n, usmin, usmax, self->state().use_early_termination);

      std::ostringstream cmd;
      
      // Use single Julia script with parameters to control early termination behavior
      cmd << "julia --project=/globalhome/tus210/HPC/quenchin_actor "
          << "/globalhome/tus210/HPC/quenchin_actor/experimental/jqsweep_custom_experimental.jl "
          << job.Ufile   << " "
          << job.Pfile   << " "
          << job.ufile   << " "
          << job.pfile   << " "
          << job.n       << " "    
          << usmax       << " "
          << job.gg      << " "
          << job.xs      << " "
          << job.theta   << " "
          << (self->state().use_early_termination ? "true" : "false") << " "  // Early termination flag
          << (self->state().use_early_termination ? "sundials" : "normal");   // Method selection

      FILE* pipe = popen(cmd.str().c_str(), "r");
      double result = 0.0;
      if (!pipe) {
        self->println("WORKER ERROR: Failed to exec Julia - popen returned NULL");
      } else {
        char buf[128];
        std::string out;
        while (fgets(buf, sizeof(buf), pipe)) out += buf;
        
        int rc = pclose(pipe);
        
        if (rc == 0) {
          try {
            result = std::stod(out);
          } catch (...) {
            self->println("WORKER ERROR: Parse error: \"{}\"", out);
          }
        } else {
          self->println("WORKER ERROR: Julia exit code {}", rc);
        }
      }

      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
      double runtime_sec = duration.count() / 1000.0;


      // send back "result", gg, theta, xs, result, usmin, usmax, runtime
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
      // spawn actor_number+1 workers
      for (int i = 0; i < actor_number + 1; ++i) {
        auto worker = self->spawn(worker_actor, self->state().manager, self->state().use_early_termination);
        // tell manager about each new worker
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
  std::unordered_map<std::tuple<int,float>, std::vector<result_entry>> uq_cache;
  std::unordered_set<std::tuple<int, float, float>> scheduled_jobs;
  std::unordered_map<std::tuple<int, float, float>, int> refinement_depth;
  
  // Spatial hash manager for efficient neighbor queries
  spatial_hash_manager spatial_mgr;

  bool use_bracket;
  bool use_early_termination;
  int  actor_number = 31;
  int  active_jobs = 0;

  // Coarse grid parameters (for adaptive refinement)
  float theta_step = 20.0f;  // Updated to match server initialization
  float xs_step;  // Will be calculated as (xs_max - xs_min) / (NDNS - 1)
  float max_theta = 40.0f;
  float xs_min = 0.1f;
  float xs_max = 700.0f;  // Maximum xs value for grid bounds
  int NDNS = 10; 

  // Grid index lookup optimization - O(1) instead of O(n) std::find_if
  std::unordered_map<float, int> theta_to_index;
  std::unordered_map<float, int> xs_to_index;

  // Worker tracking
  std::unordered_set<actor> active_workers;
  std::unordered_set<actor> idle_workers;
  
  int total_jobs_created = 0;
  int total_jobs_completed = 0;

  
  std::vector<int> gg_levels = {3, 2, 1, 0};  
  int global_refinement_level = 0;                
  std::unordered_map<int, bool> gg_level_processing;        // Whether a gg level is currently processing jobs at current refinement level
  std::unordered_map<int, bool> gg_completely_done;         // Whether gg level is completely finished (no more refinement possible)
  
  // Job tracking per gg and refinement level
  std::unordered_map<std::pair<int, int>, int> gg_level_jobs_created;    // Jobs created per (gg, refinement_level)
  std::unordered_map<std::pair<int, int>, int> gg_level_jobs_completed;  // Jobs completed per (gg, refinement_level)

  const int max_refinement_level = 4;  // Maximum refinement depth allowed
  std::chrono::high_resolution_clock::time_point start_time;
};

// Forward declarations
void bracket_initialize_from_cache(job_t& job, const std::unordered_map<std::tuple<int, float>, std::vector<result_entry>>& uq_cache, float L);

// Helper function to determine if bracket updates should be applied
bool should_apply_bracket_update(stateful_actor<server_state>* self, const job_t& job) {
  if (!self->state().use_bracket) return false;
  
  int first_gg = self->state().gg_levels[0];
  
  // Check if this is the first gg level AND a coarse grid job (depth 0)
  auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
  auto depth_it = self->state().refinement_depth.find(job_key);
  int depth = (depth_it != self->state().refinement_depth.end()) ? depth_it->second : 0;
  
  // Skip bracket updates only for coarse grid jobs (depth=0) of first gg
  if (job.gg == first_gg && depth == 0) {
    return false; // Don't apply bracket updates
  }
  
  return true; // Apply bracket updates for all other cases
}
bool is_quenching_boundary(stateful_actor<server_state>* self, const std::tuple<int, float, float>& job_key, const std::tuple<int, float, float>& neighbor_key);
void trigger_boundary_refinement(stateful_actor<server_state>* self, const std::tuple<int, float, float>& job_key, const std::tuple<int, float, float>& neighbor_key);
void schedule_gg_level_jobs(stateful_actor<server_state>* self, int gg);
void check_gg_completion_and_refine(stateful_actor<server_state>* self);
void create_all_coarse_grid_neighbors(stateful_actor<server_state>* self);
void wake_idle_workers(stateful_actor<server_state>* self);

// Grid index lookup optimization functions
void initialize_grid_lookups(stateful_actor<server_state>* self) {
  // Initialize theta lookup
  self->state().theta_to_index.clear();
  int theta_idx = 0;
  for (float theta = -self->state().max_theta; theta <= self->state().max_theta; theta += self->state().theta_step) {
    self->state().theta_to_index[theta] = theta_idx++;
  }
  
  // Initialize xs lookup  
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

// Function to create neighbors for coarse grid jobs only
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
  self->println("Scheduled job: gg={}, theta={}, xs={:.3f}, depth={}, type={}, total_created={}", 
                new_job.gg, new_job.theta, new_job.xs, depth, 
                is_refinement ? "refinement" : "coarse", self->state().total_jobs_created);
}
void check_gg_completion_and_refine(stateful_actor<server_state>* self) {
  
  bool all_gg_completed_current_level = true;
  int current_level = self->state().global_refinement_level;
  
  for (int gg : self->state().gg_levels) {
    if (self->state().gg_completely_done[gg]) continue;
    
    auto level_key = std::make_pair(gg, current_level);
    int level_created = self->state().gg_level_jobs_created[level_key];
    int level_completed = self->state().gg_level_jobs_completed[level_key];
    
    // DEBUG: Print status for each gg level
    self->println("=== GG {} LEVEL {} STATUS: jobs={}/{}, processing={}, done={} ===", 
                  gg, current_level, level_completed, level_created,
                  self->state().gg_level_processing[gg], self->state().gg_completely_done[gg]);
    
    // Check if this gg level has incomplete jobs at current level
    if (level_completed < level_created || level_created == 0) {
      all_gg_completed_current_level = false;
    }
  }
  
  // If ALL gg levels completed current level, advance ALL to next level
  if (all_gg_completed_current_level) {
    self->println("=== ALL GG LEVELS COMPLETED LEVEL {} - ADVANCING ALL TO LEVEL {} ===", 
                  current_level, current_level + 1);
    
    // Check if we can advance to next refinement level
    int next_level = current_level + 1;
    if (next_level > self->state().max_refinement_level) {
      // All levels reached maximum refinement
      self->println("=== ALL GG LEVELS REACHED MAXIMUM REFINEMENT LEVEL {} ===", self->state().max_refinement_level);
      
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
          self->state().gg_level_processing[gg] = true;
          
          std::set<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> boundary_pairs;
          std::set<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> processed_pairs;
          
          for (auto& [job_key, job_node] : self->state().job_graph) {
            auto [job_gg, job_theta, job_xs] = job_key;
            if (job_gg != gg || job_node.status != server_state::job_node::done) continue;
            
            auto neighbors_copy = job_node.neighbors;
            
            for (auto neighbor_key : neighbors_copy) {
              auto& neighbor_node = self->state().job_graph[neighbor_key];
              auto [neighbor_gg, neighbor_theta, neighbor_xs] = neighbor_key;
              
              if (neighbor_gg == gg && neighbor_node.status == server_state::job_node::done) {
                
                int depth1 = self->state().refinement_depth[job_key];
                int depth2 = self->state().refinement_depth[neighbor_key];
                int max_parent_depth = std::max(depth1, depth2);
                
                if (max_parent_depth >= self->state().max_refinement_level) {
                  continue; 
                }
                
                if (is_quenching_boundary(self, job_key, neighbor_key)) {
                  auto pair1 = std::make_pair(job_key, neighbor_key);
                  auto pair2 = std::make_pair(neighbor_key, job_key);
                  
                  if (boundary_pairs.count(pair1) == 0 && boundary_pairs.count(pair2) == 0) {
                    boundary_pairs.insert(pair1);
                  }
                }
              }
            }
          }
          
          self->println("Found {} boundary pairs for gg {} level {} → {}", 
                        boundary_pairs.size(), gg, current_level, next_level);
          
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
              
              int jobs_before = self->state().total_jobs_created;
              
              trigger_boundary_refinement(self, job_key, neighbor_key);
              
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
            self->println("=== GG {} COMPLETELY FINISHED - no boundary refinement possible ===", gg);
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
    self->println("=== ALL GG LEVELS COMPLETED ===");
    
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
  
  // Always try to wake idle workers at the end 
  wake_idle_workers(self);
}

// Helper function to wake up idle workers
void wake_idle_workers(stateful_actor<server_state>* self) {
  auto idle_it = self->state().idle_workers.begin();
  while (idle_it != self->state().idle_workers.end() && !self->state().pending_jobs.empty()) {
    auto idle_worker = *idle_it;
    idle_it = self->state().idle_workers.erase(idle_it);
    
    job_t job = self->state().pending_jobs.front();
    self->state().pending_jobs.pop_front();
    ++self->state().active_jobs;
    
    if (should_apply_bracket_update(self, job))
      bracket_initialize_from_cache(job, self->state().uq_cache, L);
    
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
  const float xs_max = self->state().xs_max;  // Use server state value
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
  
  self->println("Scheduled {} coarse grid jobs for gg level {}", jobs_for_gg, gg);
}


void bracket_initialize_from_cache(job_t& job,
  const std::unordered_map<std::tuple<int, float>, std::vector<result_entry>>& uq_cache,
  float L)
{
  // Only process quenching results for bracket updates
  bool bracket_updated = false;

  // Strategy 1: Lower gg with same theta and exact xs match
  if (!bracket_updated) {
    for (int g2 = job.gg + 1; g2 <= 3; ++g2) {
      auto it3 = uq_cache.find({g2, job.theta});
      if (it3 != uq_cache.end()) {
        for (const auto& e : it3->second) {
          if (e.uq < 0.0 && std::abs(e.xs - job.xs) < 1e-4f) {
            double adjusted_uq = e.uq * 1.5; // Increase by 50%
            
            // Apply bisection-friendly bracket adjustment: u_min = 1000/(2^n) where n = floor(log2(1000/|adjusted_uq|))
            if (adjusted_uq < 0.0) {
              double abs_adjusted_uq = std::abs(adjusted_uq);
              if (abs_adjusted_uq > 0.0) {
                job.n = static_cast<int>(std::floor(std::log2(1000.0 / abs_adjusted_uq)));
              } else {
                job.n = 0; // Fallback: usmin = -1000
              }
            } else {
              job.n = 0; // Default for non-quenching results: usmin = -1000
            }
            bracket_updated = true;
            break; // Use first match
          }
        }
        if (bracket_updated) break;
      }
    }
  }

  // Strategy 2: Exact (gg, theta) match with similar xs
  if (!bracket_updated) {
    auto it1 = uq_cache.find({job.gg, job.theta});
    if (it1 != uq_cache.end()) {
      for (const auto& e : it1->second) {
        if (e.uq < 0.0 && std::abs(e.xs - job.xs) < 0.15f * job.xs) {
          double dA = std::max(job.xs / e.xs, e.xs / job.xs);
          double adjusted_uq = e.uq * dA * 1.3; // Increase by 30% (with dA adjustment)
          
          if (adjusted_uq < 0.0) {
            double abs_adjusted_uq = std::abs(adjusted_uq);
            if (abs_adjusted_uq > 0.0) {
              job.n = static_cast<int>(std::floor(std::log2(1000.0 / abs_adjusted_uq)));
            } else {
              job.n = 0; // Fallback: usmin = -1000
            }
          } else {
            job.n = 0; // Default for non-quenching results: usmin = -1000
          }
          bracket_updated = true;
          break; // Use first match
        }
      }
    }
  }

  // Strategy 3: Same gg, different theta with exact xs match (for large xs)
  if (!bracket_updated && job.xs > (L / 10.f)) {
    for (const auto& [key, entries] : uq_cache) {
      int g_cached = std::get<0>(key);
      int t_cached = std::get<1>(key);
      if (g_cached != job.gg || t_cached == job.theta) continue;

      for (const auto& e : entries) {
        if (e.uq < 0.0 && std::abs(e.xs - job.xs) < 1e-4f) {
          double adjusted_uq = std::min(e.uq * 1.3, 0.0); // Increase by 30%
          
          // Apply bisection-friendly bracket adjustment: u_min = 1000/(2^n) where n = floor(log2(1000/|adjusted_uq|))
          if (adjusted_uq < 0.0) {
            double abs_adjusted_uq = std::abs(adjusted_uq);
            if (abs_adjusted_uq > 0.0) {
              job.n = static_cast<int>(std::floor(std::log2(1000.0 / abs_adjusted_uq)));
            } else {
              job.n = 0; // Fallback: usmin = -1000
            }
          } else {
            job.n = 0; // Default for non-quenching results: usmin = -1000
          }
          bracket_updated = true;
          break; // Use first match
        }
      }
      if (bracket_updated) break;
    }
  }
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
    self->println("*** TRUE BOUNDARY DETECTED: ({},{},{}) uq={} vs ({},{},{}) uq={} ***", 
                  std::get<0>(job_key), std::get<1>(job_key), std::get<2>(job_key), uq1,
                  std::get<0>(neighbor_key), std::get<1>(neighbor_key), std::get<2>(neighbor_key), uq2);
    is_boundary = true;
  }
  
  // For both quenching: look for significant magnitude differences
  if (!is_boundary && job_quenches && neighbor_quenches) {
    // Dynamic threshold based on gg value - different gg values have different uq ranges
    int gg = std::get<0>(job_key);  // Both jobs should have same gg value
    double threshold = 1.0;  // Default threshold
    
    // Set gg-dependent thresholds based on observed uq ranges
    if (gg == 0) {
      threshold = 0.1;  // Small threshold for gg=0
    } else if (gg == 1) {
      threshold = 0.2;  // Medium threshold for gg=1
    } else if (gg == 2) {
      threshold = 0.4;   // Larger threshold for gg=2
    } else if (gg == 3) {
      threshold = 1.0;   // Largest threshold for gg=3
    }
    
    if (std::abs(uq1 - uq2) > threshold) {
      self->println("*** QUENCH MAGNITUDE BOUNDARY: ({},{},{}) uq={} vs ({},{},{}) uq={} diff={} (threshold={}) ***", 
                    std::get<0>(job_key), std::get<1>(job_key), std::get<2>(job_key), uq1,
                    std::get<0>(neighbor_key), std::get<1>(neighbor_key), std::get<2>(neighbor_key), uq2,
                    std::abs(uq1 - uq2), threshold);
      is_boundary = true;
    }
  }
  
  // OPTIMIZATION 1: Remove this neighbor pair from both jobs' neighbor lists
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
  int max_depth = self->state().max_refinement_level; 
  auto [gg1, theta1, xs1] = job_key;
  auto [gg2, theta2, xs2] = neighbor_key;

  // Only refine within the same gg value
  if (gg1 != gg2) {
    return;  
  }

  int depth1 = self->state().refinement_depth[job_key];
  int depth2 = self->state().refinement_depth[neighbor_key];
  

  int max_parent_depth = std::max(depth1, depth2);
  if (max_parent_depth >= max_depth) {
    return;
  }

  // Calculate midpoint between the two boundary jobs
  float center_theta = (theta1 + theta2) / 2.0f;
  float center_xs = (xs1 + xs2) / 2.0f;


  // Check bounds using server state values
  if (center_theta < -self->state().max_theta || center_theta > self->state().max_theta ||
      center_xs < self->state().xs_min || center_xs > self->state().xs_max) {
    return;  
  }
  
  auto new_job_key = std::make_tuple(gg1, center_theta, center_xs);
  if (self->state().job_graph.count(new_job_key) > 0) {
    return;  
  }

  // Calculate the refinement depth for new job (increment from max parent depth)
  int new_depth = max_parent_depth + 1;
  
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
  
  schedule_job(self, new_job, true);  

  // Track refinement jobs per gg level (use global refinement level)
  int current_level = self->state().global_refinement_level;
  auto level_key = std::make_pair(gg1, current_level);
  self->state().gg_level_jobs_created[level_key]++;

  // Use depth-adaptive radius that gradually decreases but maintains connectivity
  // This balances refinement precision with neighbor connectivity
  float base_theta_spacing = self->state().theta_step;
  float base_xs_spacing = self->state().xs_step;

  float depth_factor = 1.0f / (1.0f + 0.1f * new_depth);  
  float adaptive_radius_factor = 0.6f * depth_factor;  // Base 60% scaled by depth
  float radius_theta = base_theta_spacing * adaptive_radius_factor;
  float radius_xs = base_xs_spacing * adaptive_radius_factor;
  
  // Add a small tolerance factor to ensure we catch boundary neighbors
  // This is especially important for floating-point precision issues
  float tolerance_factor = 1.2f;  // 20% extra radius for safety
  radius_theta *= tolerance_factor;
  radius_xs *= tolerance_factor;
  
  auto around = self->state().spatial_mgr.find_neighbors(new_job_key, radius_theta, radius_xs);
  
  int connections_made = 0;
  for (auto& existing_key : around) {
    auto [existing_gg, existing_theta, existing_xs] = existing_key;
    
    float d_theta = std::abs(existing_theta - center_theta);
    float d_xs = std::abs(existing_xs - center_xs);

    if (d_theta <= radius_theta && d_xs <= radius_xs) {
      self->state().job_graph[new_job_key].neighbors.insert(existing_key);
      self->state().job_graph[existing_key].neighbors.insert(new_job_key);
      connections_made++;
    }
  }
  
  self->println("REFINEMENT: Created job (gg={} global_level {} job #{}): gg={}, theta={}, xs={} (depth={}) with {} neighbors", 
                gg1, current_level, self->state().gg_level_jobs_created[level_key], gg1, center_theta, center_xs, new_depth, connections_made);
}

// Function to create complete coarse grid neighbor map after all jobs are registered
void create_all_coarse_grid_neighbors(stateful_actor<server_state>* self) {
  
  self->println("Creating coarse grid neighbor map for {} total jobs...", self->state().job_graph.size());
  
  int neighbors_created = 0;
  
  for (auto& [job_key, job_node] : self->state().job_graph) {    
    if (self->state().refinement_depth[job_key] == 0) {
      int initial_neighbors = job_node.neighbors.size();
      create_coarse_grid_neighbors(self, job_key);
      int final_neighbors = job_node.neighbors.size();
      neighbors_created += (final_neighbors - initial_neighbors);
      
      // Debug: Print neighbor count for each coarse job
      auto [gg, theta, xs] = job_key;
      self->println("Coarse job gg={}, theta={}, xs={:.1f} has {} neighbors", 
                    gg, theta, xs, final_neighbors);
    }
  }
  
  self->println("Coarse grid neighbor map complete. Total new neighbor connections: {}", neighbors_created);
}

behavior server(stateful_actor<server_state>* self,
                uint16_t port,
                bool enable_bracket,
                bool enable_early_termination) {
  self->state().use_bracket = enable_bracket;
  self->state().use_early_termination = enable_early_termination;
  self->state().start_time  = std::chrono::high_resolution_clock::now();

  if (!self->system().middleman().publish(self, port)) {
    self->println("Failed to publish server on port {}", port);
    self->quit();
  }
  self->println("BASELINE Server up on port {}, bracket-enabled={}, early_termination={}", 
                port, enable_bracket, enable_early_termination);

  // Calculate xs_step based on coarse grid parameters (only value that needs computation)
  self->state().xs_step = (self->state().xs_max - self->state().xs_min) / (self->state().NDNS - 1);

  self->state().spatial_mgr.initialize_cell_sizes(self->state().theta_step, self->state().xs_step);

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
    

    self->println("Scheduling coarse grid jobs for gg level {}...", gg);
    schedule_gg_level_jobs(self, gg);
    
    auto level_0_key = std::make_pair(gg, 0);
    auto gg_coarse_jobs = self->state().gg_level_jobs_created[level_0_key];
    
    self->println("GG {} initialized: global level 0 with {} coarse jobs", gg, gg_coarse_jobs);
  }

  // Create complete coarse grid neighbor map for ALL scheduled jobs
  create_all_coarse_grid_neighbors(self);

  self->println("=== INITIALIZATION COMPLETE: {} total jobs from {} gg levels ===", 
                self->state().pending_jobs.size(), self->state().gg_levels.size());

  for (int i = 0; i < self->state().actor_number; ++i) {
    auto worker = self->spawn(worker_actor, self, self->state().use_early_termination);
    self->state().active_workers.insert(worker);
    
    if (!self->state().pending_jobs.empty()) {
      job_t job = self->state().pending_jobs.front();
      self->state().pending_jobs.pop_front();

      if (should_apply_bracket_update(self, job))
        bracket_initialize_from_cache(job, self->state().uq_cache, L);

      auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
      self->state().job_graph[job_key].status = server_state::job_node::running;
      self->anon_send(worker, std::move(job));
      ++self->state().active_jobs;
    } else {
      // No jobs available, keep worker idle
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
          bracket_initialize_from_cache(job, self->state().uq_cache, L);
          
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
        
        self->println("*** JOB COMPLETED for gg={} global_level={}: {}/{} jobs done ***", 
                      gg, tracking_level, 
                      self->state().gg_level_jobs_completed[level_key], 
                      self->state().gg_level_jobs_created[level_key]);

        self->state().job_graph[job_key].status = server_state::job_node::done;
        self->state().job_graph[job_key].uq = uq;

        if (uq < 0.0)
            self->state().uq_cache[{gg, theta}].push_back({xs, uq});

        int neighbor_count = self->state().job_graph[job_key].neighbors.size();
        
        self->println("Job completed: gg={}, theta={}, xs={}, result={}, usmin={}, usmax={}, runtime={}s, depth={}, neighbors={}", 
                      gg, theta, xs, uq, usmin, usmax, runtime_sec, depth, neighbor_count);
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

            if (should_apply_bracket_update(self, job))
                bracket_initialize_from_cache(job, self->state().uq_cache, L);

            auto next_job_key = std::make_tuple(job.gg, job.theta, job.xs);
            self->state().job_graph[next_job_key].status = server_state::job_node::running;
            self->anon_send(worker, std::move(job));
        } else {
            // No pending jobs available - worker becomes idle
            self->state().idle_workers.insert(worker);
            
            // DEBUG: Print when worker becomes idle and check if all workers are idle
            int total_workers = self->state().active_workers.size() + self->state().idle_workers.size();
            self->println("Worker became idle. Idle workers: {}/{}, Pending jobs: {}", 
                          self->state().idle_workers.size(), total_workers, self->state().pending_jobs.size());
        }
        
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
                    check_gg_completion_and_refine(self);
                }
        }

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
        if (self->state().active_workers.empty() && self->state().idle_workers.empty()) {
          self->quit();
        }
      }
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