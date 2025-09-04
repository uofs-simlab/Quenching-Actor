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

namespace std {
  template<>
  struct hash<std::tuple<int, float, float>> {
    std::size_t operator()(const std::tuple<int, float, float>& t) const {
      auto h1 = std::hash<int>{}(std::get<0>(t));
      auto h2 = std::hash<float>{}(std::get<1>(t));
      auto h3 = std::hash<float>{}(std::get<2>(t));
      
      // Better hash combining using FNV-like approach
      std::size_t seed = 0x9e3779b9;  // Golden ratio constant
      seed ^= h1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };
  
  template<>
  struct hash<std::tuple<int, float>> {
    std::size_t operator()(const std::tuple<int, float>& t) const {
      auto h1 = std::hash<int>{}(std::get<0>(t));
      auto h2 = std::hash<float>{}(std::get<1>(t));
      
      std::size_t seed = 0x9e3779b9;
      seed ^= h1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };
  
  template<>
  struct hash<std::tuple<int, int, int>> {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
      auto h1 = std::hash<int>{}(std::get<0>(t));
      auto h2 = std::hash<int>{}(std::get<1>(t));
      auto h3 = std::hash<int>{}(std::get<2>(t));
      
      std::size_t seed = 0x9e3779b9;
      seed ^= h1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };
  
  template<>
  struct hash<std::pair<int, int>> {
    std::size_t operator()(const std::pair<int, int>& p) const {
      auto h1 = std::hash<int>{}(p.first);
      auto h2 = std::hash<int>{}(p.second);
      
      std::size_t seed = 0x9e3779b9;
      seed ^= h1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };
  
  template<>
  struct hash<std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>> {
    std::size_t operator()(const std::pair<std::tuple<int, float, float>, std::tuple<int, float, float>>& p) const {
      auto h1 = std::hash<std::tuple<int, float, float>>{}(p.first);
      auto h2 = std::hash<std::tuple<int, float, float>>{}(p.second);
      
      std::size_t seed = 0x9e3779b9;
      seed ^= h1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };
} 

struct job_t {
  int      gg;
  float    theta;
  float    xs;
  int      n;           // Power for bracket: usmin = -1000/(2^n), usmax = 0
  std::string Ufile;
  std::string Pfile;
  std::string ufile;
  std::string pfile;
  // std::string savedir;

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
      
      cmd << "julia --project=/globalhome/tus210/HPC/quenchin_actor "
          << "/globalhome/tus210/HPC/quenchin_actor/experimental/jqsweep_custom_experimental.jl "
          << job.Ufile   << " "
          << job.Pfile   << " "
          << job.ufile   << " "
          << job.pfile   << " "
          << job.n       << " "    // Pass n instead of usmin
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
        
        // self->println("WORKER DEBUG: Julia output: '{}'", out);
        
        int rc = pclose(pipe);
        // self->println("WORKER DEBUG: Julia exit code: {}", rc);
        
        if (rc == 0) {
          try {
            result = std::stod(out);
            // self->println("-> result = {}", result);
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

      // send back "result", gg, theta, xs, result, runtime
      anon_mail("result", job.gg, job.theta, job.xs, result, usmin, usmax,
         runtime_sec, self)
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

// Optimized spatial hash manager for efficient neighbor queries
struct spatial_hash_manager {
    float theta_cell_size = 2.0f;  // Will be set based on coarse grid parameters
    float xs_cell_size = 10.0f;    // Will be set based on coarse grid parameters
    
    // Initialize cell sizes based on coarse grid spacing
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
  float theta_step = 10.0f;
  float xs_step = 50.0f;  // Will be calculated from xs_max, xs_min, NDNS
  float max_theta = 40.0f;
  float xs_min = 0.1f;
  int NDNS = 10;  // Match baseline: 10 xs values, not 15

  // Pre-computed grid index lookups to avoid O(n) std::find_if operations
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

  std::chrono::high_resolution_clock::time_point start_time;
};



// Initialize grid index lookups for O(1) grid index queries
void initialize_grid_lookups(stateful_actor<server_state>* self) {
  self->state().theta_to_index.clear();
  self->state().xs_to_index.clear();
  self->state().coarse_theta_values.clear();
  self->state().coarse_xs_values.clear();
  
  float theta_step = self->state().theta_step;
  float xs_step = self->state().xs_step;
  float max_theta = self->state().max_theta;
  float xs_min = self->state().xs_min;
  int NDNS = self->state().NDNS;
  
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
  return -1; // Not found on coarse grid
}

int get_xs_index(stateful_actor<server_state>* self, float xs) {
  // Find closest xs value in coarse grid  
  for (const auto& [grid_xs, index] : self->state().xs_to_index) {
    if (std::abs(grid_xs - xs) < 1.0f) {
      return index;
    }
  }
  return -1; // Not found on coarse grid
}

// Function to create neighbors for coarse grid jobs only
void create_coarse_grid_neighbors(stateful_actor<server_state>* self, 
                                 const std::tuple<int, float, float>& job_key) {
  
  auto [gg, thetaC, xsC] = job_key;
  
  // Get grid indices with O(1) lookup instead of O(n) std::find_if
  int theta_idx = get_theta_index(self, thetaC);
  int xs_idx = get_xs_index(self, xsC);
  
  if (theta_idx == -1 || xs_idx == -1) {
    return; // Not a coarse grid job
  }

  // Only connect to other coarse grid jobs
  for (auto& [keyE, nodeE] : self->state().job_graph) {
      if (keyE == job_key) continue; // skip self
      if (std::get<0>(keyE) != gg) continue; // only same gg
      
      // Only connect to other coarse jobs (refinement_depth == 0)
      if (self->state().refinement_depth[keyE] != 0) continue;

      auto [_, thetaE, xsE] = keyE;

      // Get grid indices with O(1) lookup instead of O(n) std::find_if
      int theta_idx_E = get_theta_index(self, thetaE);
      int xs_idx_E = get_xs_index(self, xsE);
      
      if (theta_idx_E == -1 || xs_idx_E == -1) {
        continue; // Other job not on coarse grid
      }
      
      // EXPLICIT 8-NEIGHBOR CONNECTIVITY (includes all diagonals)
      // Calculate grid distance for both dimensions
      int delta_theta = std::abs(theta_idx - theta_idx_E);
      int delta_xs = std::abs(xs_idx - xs_idx_E);
      
      // Connect to all 8 surrounding neighbors:
      // Horizontal: (±1,0), Vertical: (0,±1), Diagonal: (±1,±1)
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
}


void bracket_initialize_from_cache(job_t& job,
  const std::unordered_map<std::tuple<int, float>, std::vector<result_entry>>& uq_cache,
  float L,
  stateful_actor<server_state>* self = nullptr)
{
  // Only process quenching results for bracket updates
  bool bracket_updated = false;

  // STRATEGY 1: Check neighbors that are done and update usmin to their uq with trust margin
  if (!bracket_updated && self != nullptr) {
    try {
      auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
      
      // Check if this job exists in job_graph and has neighbors
      auto job_graph_it = self->state().job_graph.find(job_key);
      if (job_graph_it != self->state().job_graph.end() && !job_graph_it->second.neighbors.empty()) {
        auto& node = job_graph_it->second;
        
        // Make a copy of neighbors to avoid iterator invalidation
        auto neighbors_copy = node.neighbors;
        
        for (auto neighbor_key : neighbors_copy) {
          auto neighbor_it = self->state().job_graph.find(neighbor_key);
          if (neighbor_it != self->state().job_graph.end() && 
              neighbor_it->second.status == server_state::job_node::done && 
              neighbor_it->second.uq < 0.0) {
            
            // Found a done neighbor with quenching result - use its uq with trust margin
            double adjusted_uq = neighbor_it->second.uq * 1.5; // Trust margin of 50%
            
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
            break; // Use first neighbor match
          }
        }
      }
    } catch (...) {
      // Silently skip Strategy 1 if there are any issues
    }
  }

  // STRATEGY 2: Lower gg with same theta and exact xs match (from cache)
  if (!bracket_updated) {
    for (int g2 = job.gg + 1; g2 <= 3; ++g2) {
      auto it3 = uq_cache.find({g2, job.theta});
      if (it3 != uq_cache.end()) {
        for (const auto& e : it3->second) {
          if (e.uq < 0.0 && std::abs(e.xs - job.xs) < 1e-4f) {
            double adjusted_uq = e.uq * 1.5; // Trust margin of 50%
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
}

bool is_quenching_boundary(stateful_actor<server_state>* self,
                           const std::tuple<int, float, float>& job_key,
                           const std::tuple<int, float, float>& neighbor_key) {
  
  auto& job_node = self->state().job_graph[job_key];
  auto& neighbor_node = self->state().job_graph[neighbor_key];
  
  double uq1 = job_node.uq;
  double uq2 = neighbor_node.uq;
  
  // Primary boundary criterion: one quenches, one doesn't
  bool job_quenches = (uq1 < 0.0);
  bool neighbor_quenches = (uq2 < 0.0);
  bool is_boundary = false;
  
  if (job_quenches != neighbor_quenches) {
    is_boundary = true;
  }
  
  // For both quenching: look for significant magnitude differences
  if (!is_boundary && job_quenches && neighbor_quenches) {
    // Dynamic threshold based on gg value - different gg values have different uq ranges
    int gg = std::get<0>(job_key);  // Both jobs should have same gg value
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

// Helper function to get the maximum refinement depth
int get_max_refinement_depth() {
  return 2; // Centralized max depth configuration
}

void trigger_boundary_refinement(stateful_actor<server_state>* self,
                                const std::tuple<int, float, float>& job_key,
                                const std::tuple<int, float, float>& neighbor_key) {
  static const int max_depth = get_max_refinement_depth(); // Use centralized max depth
  auto [gg1, theta1, xs1] = job_key;
  auto [gg2, theta2, xs2] = neighbor_key;

  // Only refine within the same gg value
  if (gg1 != gg2) {
    return;
  }

  // Check refinement depth for both jobs
  int depth1 = self->state().refinement_depth[job_key];
  int depth2 = self->state().refinement_depth[neighbor_key];
  
  // Check if we would create jobs beyond max_depth
  int max_parent_depth = std::max(depth1, depth2);
  if (max_parent_depth >= max_depth) {
    return;
  }

  // Calculate midpoint between the two boundary jobs
  float center_theta = (theta1 + theta2) / 2.0f;
  float center_xs = (xs1 + xs2) / 2.0f;

  // Check bounds
  if (center_theta < -40.0f || center_theta > 40.0f ||
      center_xs < 0.1f || center_xs > 700.0f) {
    return;
  }
  
  auto new_job_key = std::make_tuple(gg1, center_theta, center_xs);
  if (self->state().job_graph.count(new_job_key) > 0) {
    return;
  }

  // Calculate the refinement depth for new job (increment from max parent depth)
  int new_depth = max_parent_depth + 1;
  
  // Calculate base spacing early for use in minimum distance calculations
  float base_theta_spacing = self->state().theta_step;  // 10.0f
  float base_xs_spacing = self->state().xs_step;        // calculated from coarse grid
  
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
  
  // Schedule the job to add it to spatial hash and job graph (add to END of queue)
  schedule_job(self, new_job, true);  // true = is_refinement

  // OPTIMIZATION 2: ADAPTIVE NEIGHBOR RADIUS for better distribution
  // Use depth-adaptive radius that gradually decreases but maintains connectivity
  // This balances refinement precision with neighbor connectivity
  
  // ADAPTIVE RADIUS: Slightly decrease with depth but maintain good connectivity
  // This balances connectivity with performance at deeper levels
  float depth_factor = 1.0f / (1.0f + 0.1f * new_depth);  // Gradual decrease: 1.0, 0.91, 0.83, 0.77...
  float fixed_radius_factor = 0.6f * depth_factor;  // Start at 60%, gradually decrease
  float radius_theta = base_theta_spacing * fixed_radius_factor;
  float radius_xs = base_xs_spacing * fixed_radius_factor;
  
  // Add a small tolerance factor to ensure we catch boundary neighbors
  // This is especially important for floating-point precision issues
  float tolerance_factor = 1.2f;  // 20% extra radius for safety
  radius_theta *= tolerance_factor;
  radius_xs *= tolerance_factor;
  
  // Optional: Add minimum radius to ensure connectivity
  radius_theta = std::max(radius_theta, base_theta_spacing * 0.3f);
  radius_xs = std::max(radius_xs, base_xs_spacing * 0.3f);

  auto around = self->state().spatial_mgr.find_neighbors(new_job_key, radius_theta, radius_xs);
  
  // Connect to all jobs within the fixed radius
  // This naturally includes horizontal, vertical, and diagonal connections
  int connections_made = 0;
  for (auto& existing_key : around) {
    auto [existing_gg, existing_theta, existing_xs] = existing_key;
    
    // Note: find_neighbors already filters by same gg, so no need to check again

    float d_theta = std::abs(existing_theta - center_theta);
    float d_xs = std::abs(existing_xs - center_xs);

    // This creates natural 8-directional connectivity (horizontal, vertical, diagonal)
    if (d_theta <= radius_theta && d_xs <= radius_xs) {
      self->state().job_graph[new_job_key].neighbors.insert(existing_key);
      self->state().job_graph[existing_key].neighbors.insert(new_job_key);
      connections_made++;
    }
  }
  
  self->println("REFINEMENT: Created job: gg={}, theta={}, xs={} (depth={}) with {} neighbors", 
                gg1, center_theta, center_xs, new_depth, connections_made);
}

// Function to create complete coarse grid neighbor map after all jobs are registered
void create_all_coarse_grid_neighbors(stateful_actor<server_state>* self) {
  
  // Create all neighbor relationships for coarse grid jobs
  for (auto& [job_key, job_node] : self->state().job_graph) {    
    // Only process coarse grid jobs
    if (self->state().refinement_depth[job_key] == 0) {
      create_coarse_grid_neighbors(self, job_key);
    }
  }
  
  // Count total neighbors for verification
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
  self->state().use_bracket = enable_bracket;
  self->state().use_early_termination = enable_early_termination;
  self->state().start_time  = std::chrono::high_resolution_clock::now();

  if (!self->system().middleman().publish(self, port)) {
    self->println("Failed to publish server on port {}", port);
    self->quit();
  }
  self->println("Server up on port {}, bracket-enabled={}, early-termination={}", 
                port, enable_bracket, enable_early_termination);

  // Use the grid parameters already defined in server_state struct
  float L = 2700.f;
  const float xs_max = 700.0f;

  // No need to redefine - use state values directly
  // self->state().theta_step is already 10.0f
  // self->state().max_theta is already 40.0f  
  // self->state().xs_min is already 0.1f
  // self->state().NDNS is already 10

  // Calculate xs_step using state values
  self->state().xs_step = (xs_max - self->state().xs_min) / (self->state().NDNS - 1);
  
  // Initialize spatial hash cell sizes based on coarse grid parameters
  self->state().spatial_mgr.initialize_cell_sizes(self->state().theta_step, self->state().xs_step);

  // OPTIMIZATION: Initialize grid index lookups for O(1) access
  initialize_grid_lookups(self);

  std::vector<float> xs_values;
  for (int i = 0; i < self->state().NDNS; ++i) {
    float xs = self->state().xs_min + i * self->state().xs_step;
    xs_values.push_back(xs);
  }

  std::vector<int> all_theta;
  for (float theta = -self->state().max_theta; theta <= self->state().max_theta; theta += self->state().theta_step)
    all_theta.push_back(theta);

  // OPTIMIZATION: Pre-calculate job count and reserve hash map capacity to avoid rehashing
  // Number of jobs = gg_count * theta_count * xs_count
  std::vector<int> test_gg_values = {3, 2, 1, 0};  // Process gg=3 first, down to gg=0
  size_t estimated_job_count = test_gg_values.size() * all_theta.size() * xs_values.size();
  size_t estimated_capacity = static_cast<size_t>(estimated_job_count * 1.5); // Extra capacity for refinement jobs
  
  // Reserve capacity in major hash maps to avoid expensive rehashing during job creation
  self->state().job_graph.reserve(estimated_capacity);
  self->state().scheduled_jobs.reserve(estimated_capacity);
  self->state().refinement_depth.reserve(estimated_capacity);
  self->state().spatial_mgr.grid.reserve(estimated_capacity / 10); // Spatial cells will be much fewer
  
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
      // Track the worker
      self->state().active_workers.insert(worker);
      
      if (!self->state().pending_jobs.empty()) {
        job_t job = self->state().pending_jobs.front();
        self->state().pending_jobs.pop_front();
        if (self->state().use_bracket)
          bracket_initialize_from_cache(job, self->state().uq_cache, L, self);
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

        // Update job_graph
        self->state().job_graph[job_key].status = server_state::job_node::done;
        self->state().job_graph[job_key].uq = uq;

        // Update cache if uq < 0
        if (uq < 0.0)
            self->state().uq_cache[{gg, theta}].push_back({xs, uq});

        int depth = 0;
        auto depth_it = self->state().refinement_depth.find(job_key);
        if (depth_it != self->state().refinement_depth.end())
          depth = depth_it->second;
        
        // Get neighbor count for this job
        int neighbor_count = self->state().job_graph[job_key].neighbors.size();

        self->println("Job completed: gg={}, theta={}, xs={}, result={}, usmin={}, usmax={}, runtime={}s, depth={}, neighbors={}", 
                      gg, theta, xs, uq, usmin, usmax, runtime_sec, depth, neighbor_count);

        // === AUTOMATIC RESULT DEDUCTION ===
        // COMMENTED OUT: Auto-deduction creates inconsistency with boundary detection
        // Jobs that are deduced (not actually run) don't trigger boundary refinement
        // This causes different refinement patterns compared to baseline
        /*
        // If result is 1.0 (not quenched), all jobs with same (gg, theta) but smaller xs are also 1.0
        if (uq >= 0.0) {  // Job did not quench (result = 1.0)
            // Check all pending jobs
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
                    int pending_depth = 0;
                    auto pending_depth_it = self->state().refinement_depth.find(pending_key);
                    if (pending_depth_it != self->state().refinement_depth.end())
                      pending_depth = pending_depth_it->second;
                    int pending_neighbor_count = self->state().job_graph[pending_key].neighbors.size();
                    
                    // Print the result like a normal job completion
                    self->println("Job completed: gg={}, theta={}, xs={}, result={}, usmin={}, usmax={}, runtime={}s, depth={}, neighbors={}", 
                                  pending_job.gg, pending_job.theta, pending_job.xs, 1.0, usmin, usmax, 0.0, pending_depth, pending_neighbor_count);

                    // Remove from scheduled jobs
                    self->state().scheduled_jobs.erase(pending_key);
                    
                    // Increment completion counter
                    ++self->state().total_jobs_completed;
                    
                    // Remove from pending queue
                    it = self->state().pending_jobs.erase(it);
                  }
                  else {
                    ++it;
                  }
              }
          }
        */

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

            if (self->state().use_bracket)
                bracket_initialize_from_cache(job, self->state().uq_cache, L, self);

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
        
        // First pass: identify all boundary pairs with pre-filtering
        for (auto neighbor_key : neighbors_copy) {
            auto& neighbor_node = self->state().job_graph[neighbor_key];

            // Only process neighbors that are already done
            if (neighbor_node.status == server_state::job_node::done) {
                
                // OPTIMIZATION: Pre-filter boundary pairs that have already reached max depth
                // This prevents calling trigger_boundary_refinement for pairs that will be rejected
                int depth1 = self->state().refinement_depth[job_key];
                int depth2 = self->state().refinement_depth[neighbor_key];
                int max_parent_depth = std::max(depth1, depth2);
                
                if (max_parent_depth >= get_max_refinement_depth()) {  // Use centralized max depth
                    // Skip - this pair has already reached maximum refinement depth
                    continue;
                }
                
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
        
        // Second pass: process all identified boundary pairs efficiently
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

        // === Check if new jobs were created and wake idle workers ===
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
                
                if (self->state().use_bracket)
                    bracket_initialize_from_cache(job, self->state().uq_cache, L, self);
                
                auto wake_job_key = std::make_tuple(job.gg, job.theta, job.xs);
                self->state().job_graph[wake_job_key].status = server_state::job_node::running;
                self->anon_send(idle_worker, std::move(job));
                
                // self->println("WAKE: Sent refinement job to previously idle worker {}", idle_worker.id());
            }
        }

        // === Check for completion ===
        if (self->state().pending_jobs.empty() && self->state().active_jobs == 0) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - self->state().start_time);
            self->println("=== ALL JOBS COMPLETED ===");
            self->println("Total runtime: {} seconds", duration.count());
            self->println("Jobs completed: {}/{}", self->state().total_jobs_completed, self->state().total_jobs_created);
            
            // Send quit message to all workers
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