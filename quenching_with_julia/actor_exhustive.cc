#include "config.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cmath>

// JOB & CACHE TYPES
struct job_t {
  int      gg;
  int      theta;
  float    xs;
  int      n;           // Power for bracket: usmin = -1000/(2^n), usmax = 0
  std::string Ufile;
  std::string Pfile;
  std::string ufile;
  std::string pfile;
  std::string savedir;

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
      f.field("pfile",  x.pfile),
      f.field("savedir",x.savedir)
    );
  }
};

CAF_BEGIN_TYPE_ID_BLOCK(my_project, caf::first_custom_type_id)
CAF_ADD_TYPE_ID(my_project, (job_t))
CAF_END_TYPE_ID_BLOCK(my_project)

// Forward declarations
struct running_job_info;
struct server_state;
double estimate_job_runtime(const job_t& job, const std::map<std::tuple<int, int, float>, double>& result_map);
bool should_restart_job(const running_job_info& job_info, int new_strategy, 
                       const std::map<std::tuple<int, int, float>, double>& result_map);
void check_for_potential_bracket_updates(stateful_actor<server_state>* self, float L);
void generate_comprehensive_timing_report(stateful_actor<server_state>* self, float L);


// WORKER
struct worker_state {
  actor manager_actor;
  bool use_early_termination;
};
behavior worker_actor(stateful_actor<worker_state>* self, actor manager, bool use_early_termination) {
  self->state().manager_actor = manager;
  self->state().use_early_termination = use_early_termination;
  self->println("Worker spawned");

  return {
    [=](job_t job) {
      auto job_start_time = std::chrono::high_resolution_clock::now();
      
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
        self->println("Failed to exec Julia");
      } else {
        char buf[128];
        std::string out;
        while (fgets(buf, sizeof(buf), pipe)) out += buf;
        int rc = pclose(pipe);
        if (rc == 0) {
          try {
            result = std::stod(out);
            // self->println("-> result = {}", result);
          } catch (...) {
            self->println("Parse error: \"{}\"", out);
          }
        } else {
          self->println("Julia exit code {}", rc);
        }
      }

      // Calculate job execution time
      auto job_end_time = std::chrono::high_resolution_clock::now();
      auto job_duration = std::chrono::duration_cast<std::chrono::milliseconds>(job_end_time - job_start_time);
      double job_runtime_seconds = job_duration.count() / 1000.0;

      std::ofstream logfile("/globalhome/tus210/HPC/quenchin_actor/caf_runtime_log.txt", std::ios::app);
      if (logfile.is_open()) {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        auto tm = *std::localtime(&time_t);
        
        // Compute usmin from n for logging
        double usmin = -1000.0 / std::pow(2.0, job.n);
        
        logfile << "[" << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S") << "] "
                << "gg=" << job.gg << " "
                << "theta=" << job.theta << " "
                << "xs=" << job.xs << " "
                << "n=" << job.n << " "
                << "usmin=" << usmin << " "
                << "usmax=" << 0.0 << " "
                << "runtime=" << job_runtime_seconds << "s "
                << "result=" << result << std::endl;
        logfile.close();
      }
      
      self->println("Job completed in {}s: gg={}, theta={}, xs={}, n={}, usmin={}, usmax={}, result={}", 
                    job_runtime_seconds, job.gg, job.theta, job.xs, job.n, usmin, usmax, result);

      // send back "result", gg, theta, xs, result, execution_time
      anon_mail("result", job.gg, job.theta, job.xs, result, job_runtime_seconds, self)
        .send(self->state().manager_actor);
    },

    [=](const std::string& msg) {
      if (msg == "quit") {
        self->quit();
      }
    },
  };
}

// REMOTE
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
        // tell manager about each new worker
        anon_mail(worker).send(self->state().manager);
      }
      self->quit();
    }
  };
}


/// SERVER

struct running_job_info {
  job_t job;
  actor worker;
  std::chrono::high_resolution_clock::time_point start_time;
  int current_strategy;
};

struct server_state {
  std::vector<job_t> default_jobs;     // Jobs with default bracket [-1000, 0]
  std::vector<job_t> updated_jobs;     // Jobs with updated brackets (higher priority)
  std::map<std::tuple<int, int, float>, double> result_map;
  std::map<std::tuple<int, int, float>, int> strategy_map; 
  
  // Track running jobs with their start times
  std::map<actor, running_job_info> running_jobs;

  bool use_bracket;
  bool use_early_termination;
  int  actor_number = 31;
  int  jobs_done    = 0;
  int  active_jobs  = 0;

  std::chrono::high_resolution_clock::time_point start_time;
};

// Function to update default jobs with better brackets based on a new result
void update_jobs_with_result(server_state& state, int completed_gg, int completed_theta, 
                             float completed_xs, double completed_uq, float L,
                             std::function<void(const std::string&)> logger = nullptr) {
  if (completed_uq >= 0.0) return; 
  
  auto it = state.default_jobs.begin();
  int updates_applied = 0;
  
  while (it != state.default_jobs.end()) {
    job_t& job = *it;
    bool should_update = false;
    int strategy_applied = 0;
    
    // Strategy 1: Lower gg with same theta and exact xs match
    if (!should_update && job.gg < completed_gg && job.theta == completed_theta) {
      if (std::abs(job.xs - completed_xs) < 1e-4f) {
        double adjusted_uq = completed_uq * 1.5; // Increase by 50% 
        if (adjusted_uq < 0.0) {
          double abs_adjusted_uq = std::abs(adjusted_uq);
          if (abs_adjusted_uq > 0.0) {
            int n = static_cast<int>(std::floor(std::log2(1000.0 / abs_adjusted_uq)));
            job.n = n;
          } else {
            job.n = 4; 
          }
        } else {
          job.n = 4; 
        }
        strategy_applied = 3;
        should_update = true;
        
        if (logger) {
          double usmin_calc = -1000.0 / std::pow(2.0, job.n);
          logger("Strategy 3 update: job_gg=" + std::to_string(job.gg) + ", theta=" + std::to_string(job.theta) + 
                 ", xs=" + std::to_string(job.xs) + " updated from completed_gg=" + std::to_string(completed_gg) + 
                 ", uq=" + std::to_string(completed_uq) + ", adjusted_uq=" + std::to_string(adjusted_uq) +
                 ", n=" + std::to_string(job.n) + ", new bracket=[" + std::to_string(usmin_calc) + 
                 ", 0.0]");
        }
      }
    }

    // Strategy 1: Exact (gg, theta) match with similar xs
    if (job.gg == completed_gg && job.theta == completed_theta) {
      if (std::abs(job.xs - completed_xs) < 0.15f * job.xs) {
        double dA = std::max(job.xs / completed_xs, completed_xs / job.xs);
        double adjusted_uq = completed_uq * dA * 1.3; // Increase by 30% (with dA adjustment) instead of multiplying by 10
        
        // Apply bisection-friendly bracket adjustment: u_min = 1000/(2^n) where n = floor(log2(1000/|adjusted_uq|))
        if (adjusted_uq < 0.0) {
          double abs_adjusted_uq = std::abs(adjusted_uq);
          if (abs_adjusted_uq > 0.0) {
            int n = static_cast<int>(std::floor(std::log2(1000.0 / abs_adjusted_uq)));
            job.n = n;
          } else {
            job.n = 4; // Fallback for edge case (equivalent to usmin = -62.5)
          }
        } else {
          job.n = 4; // Default for non-quenching results
        }
        strategy_applied = 1;
        should_update = true;
        
        if (logger) {
          double usmin_calc = -1000.0 / std::pow(2.0, job.n);
          logger("Strategy 1 update: gg=" + std::to_string(job.gg) + ", theta=" + std::to_string(job.theta) + 
                 ", xs=" + std::to_string(job.xs) + " updated from completed xs=" + std::to_string(completed_xs) + 
                 ", uq=" + std::to_string(completed_uq) + ", adjusted_uq=" + std::to_string(adjusted_uq) +
                 ", n=" + std::to_string(job.n) + ", new bracket=[" + std::to_string(usmin_calc) + 
                 ", 0.0]");
        }
      }
    }

    
    // Strategy 2: Same gg, different theta with exact xs match (for large xs)
    if (!should_update && job.gg == completed_gg && job.theta != completed_theta && job.xs > (L / 10.f)) {
      if (std::abs(job.xs - completed_xs) < 1e-4f) {
        double adjusted_uq = std::min(completed_uq * 1.3, 0.0); // Increase by 30% instead of multiplying by 10
        
        // Apply bisection-friendly bracket adjustment: u_min = 1000/(2^n) where n = floor(log2(1000/|adjusted_uq|))
        if (adjusted_uq < 0.0) {
          double abs_adjusted_uq = std::abs(adjusted_uq);
          if (abs_adjusted_uq > 0.0) {
            int n = static_cast<int>(std::floor(std::log2(1000.0 / abs_adjusted_uq)));
            job.n = n;
          } else {
            job.n = 4; // Fallback for edge case (equivalent to usmin = -62.5)
          }
        } else {
          job.n = 4; // Default for non-quenching results
        }
        strategy_applied = 2;
        should_update = true;
        
        if (logger) {
          double usmin_calc = -1000.0 / std::pow(2.0, job.n);
          logger("Strategy 2 update: gg=" + std::to_string(job.gg) + ", job_theta=" + std::to_string(job.theta) + 
                 ", xs=" + std::to_string(job.xs) + " updated from completed_theta=" + std::to_string(completed_theta) + 
                 ", uq=" + std::to_string(completed_uq) + ", adjusted_uq=" + std::to_string(adjusted_uq) +
                 ", n=" + std::to_string(job.n) + ", new bracket=[" + std::to_string(usmin_calc) + 
                 ", 0.0]");
        }
      }
    }
    
    if (should_update) {
      // Store the strategy for this job
      state.strategy_map[{job.gg, job.theta, job.xs}] = strategy_applied;
      
      // Move the updated job to the updated_jobs queue (higher priority)
      state.updated_jobs.push_back(*it);
      it = state.default_jobs.erase(it);
      updates_applied++;
    } else {
      ++it;
    }
  }
  
  if (logger && updates_applied > 0) {
    logger("Applied " + std::to_string(updates_applied) + " bracket updates. " + 
           std::to_string(state.default_jobs.size()) + " jobs remain with default brackets.");
  }
}

static std::vector<float> make_uniform_xs(float xs_min,
                                          float xs_max,
                                          int count) {
  std::vector<float> xs;
  xs.reserve(count);
  if (count == 1) {
    xs.push_back(0.5f*(xs_min + xs_max));
  } else {
    float dx = (xs_max - xs_min) / float(count - 1);
    for (int i = 0; i < count; ++i) {
      xs.push_back(xs_min + dx * i);
    }
  }
  return xs;
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
  self->println("Server up on port {}, bracket-enabled={}, early_termination={}", 
                port, enable_bracket, enable_early_termination);

  float L = 2700.f;
  int NDNS = 15;
  int max_theta = 40;
  int theta_step = 10;
  const float xs_min = 0.1f;
  const float xs_max = 700.0f;
  const float dx = L / NDNS;

  std::vector<int> all_theta;
  for (int theta = -max_theta; theta <= max_theta; theta += theta_step)
    all_theta.push_back(theta);

  for (int gg = 0; gg <= 3; ++gg) {

    std::vector<float> xs_values = make_uniform_xs(xs_min, xs_max, NDNS);


    std::string base = "/globalhome/tus210/HPC/quenchin_actor/";
    std::string waved = base + "waves/index_11";
    std::string crit = base + "waves/index_10";
    auto Ufile = waved + "/" + std::to_string(gg) + "/U";
    auto Pfile = waved + "/" + std::to_string(gg) + "/p";
    auto ufile = crit + "/" + std::to_string(gg) + "/U";
    auto pfile = crit + "/" + std::to_string(gg) + "/p";

    for (float xs : xs_values) {
      for (int theta : all_theta) {
        job_t job{
          gg, theta, xs, 0,  // n = 0 gives usmin = -1000/(2^0) = -1000, usmax = 0
          Ufile, Pfile, ufile, pfile,
          base + "dns/" + std::to_string(gg) + "/" + std::to_string(theta) + "/xs_" + std::to_string(xs) + "/"
        };
        self->state().default_jobs.push_back(job);
      }
    }
  }


  self->println("Enqueued {} jobs (all with default n=0, equivalent to bracket [-1000, 0])", self->state().default_jobs.size());

  for (int i = 0; i < self->state().actor_number; ++i) {
    auto worker = self->spawn(worker_actor, self, self->state().use_early_termination);
    if (!self->state().default_jobs.empty()) {
      job_t job = self->state().default_jobs.back();
      self->state().default_jobs.pop_back();
      
      // Track running jobs for runtime calculation
      running_job_info job_info;
      job_info.job = job;
      job_info.worker = worker;
      job_info.start_time = std::chrono::high_resolution_clock::now();
      job_info.current_strategy = 0; // Default strategy (no bracket update)
      
      self->state().running_jobs[worker] = job_info;
      
      ++self->state().active_jobs;
      self->anon_send(worker, std::move(job));
    } else {
      self->anon_send(worker, "quit");
    }
  }

  return {
    [=](actor remote, int) {
      anon_mail(self->state().actor_number).send(remote);
    },

    [=](actor worker) {
      job_t job;
      bool has_job = false;
      bool is_updated_job = false;
      
      // First priority: get a job with updated brackets
      if (!self->state().updated_jobs.empty()) {
        job = self->state().updated_jobs.back();
        self->state().updated_jobs.pop_back();
        has_job = true;
        is_updated_job = true;
      } 
      // Second priority: get a job with default brackets
      else if (!self->state().default_jobs.empty()) {
        job = self->state().default_jobs.back();
        self->state().default_jobs.pop_back();
        has_job = true;
        is_updated_job = false;
      }
      
      if (has_job) {
        running_job_info job_info;
        job_info.job = job;
        job_info.worker = worker;
        job_info.start_time = std::chrono::high_resolution_clock::now();
        
        if (is_updated_job) {
          // Get the strategy that was applied to this job
          auto strategy_it = self->state().strategy_map.find({job.gg, job.theta, job.xs});
          job_info.current_strategy = (strategy_it != self->state().strategy_map.end()) ? strategy_it->second : 0;
          double usmin_calc = -1000.0 / std::pow(2.0, job.n);
          self->println("Assigning updated job with strategy {} and bracket [{}, {}] for gg={}, theta={}, xs={}", 
                        job_info.current_strategy, usmin_calc, 0.0, job.gg, job.theta, job.xs);
        } else {
          job_info.current_strategy = 0; // Default strategy (no bracket update)
          double usmin_calc = -1000.0 / std::pow(2.0, job.n);
          self->println("Assigning default job with bracket [{}, {}] for gg={}, theta={}, xs={}", 
                        usmin_calc, 0.0, job.gg, job.theta, job.xs);
        }
        
        self->state().running_jobs[worker] = job_info;
        ++self->state().active_jobs;
        self->anon_send(worker, std::move(job));
      } else {
        self->anon_send(worker, "quit");
      }
    },

    [=](const std::string& msg, int gg, int theta, float xs, double uq, double job_runtime_seconds, actor worker) {
      if (msg != "result") return;

      
      job_t completed_job; 
      bool job_found = false;
      
      auto running_it = self->state().running_jobs.find(worker);
      if (running_it != self->state().running_jobs.end()) {
        completed_job = running_it->second.job;
        job_found = true;
        

        self->state().running_jobs.erase(worker);
      } else {
        self->println("Warning: Job info not found for worker");
      }

      std::ofstream runtime_log;
      std::string log_filename = self->state().use_bracket ? 
        "caf_exhustive_runtime_log.txt" : "caf_exhustive_without_runtime_log.txt";
      runtime_log.open(log_filename, std::ios::app);
      if (runtime_log.is_open()) {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        
        std::stringstream timestamp;
        timestamp << std::put_time(std::gmtime(&time_t), "%Y-%m-%dT%H:%M:%S");
        
        runtime_log << "[" << timestamp.str() << "] "
                   << "gg=" << gg << " "
                   << "theta=" << theta << " "
                   << "xs=" << xs << " ";
        
        if (job_found) {
          double usmin_calc = -1000.0 / std::pow(2.0, completed_job.n);
          runtime_log << "usmin=" << usmin_calc << " usmax=0.0 ";
        } else {
          runtime_log << "usmin=unknown usmax=unknown ";
        }
        
        runtime_log << "execution_time=" << job_runtime_seconds << "s "
                   << "result=" << uq << std::endl;
        runtime_log.close();
      }

      self->state().jobs_done++;
      --self->state().active_jobs;
      self->state().result_map[{gg, theta, xs}] = uq;

      if (self->state().use_bracket && uq < 0.0) {
        auto bracket_start = std::chrono::high_resolution_clock::now();
        
        auto logger = [](const std::string& msg) {
          std::ofstream bracket_log("bracket_strategy_log.txt", std::ios::app);
          if (bracket_log.is_open()) {
            auto now = std::chrono::system_clock::now();
            auto time_t = std::chrono::system_clock::to_time_t(now);
            std::stringstream timestamp;
            timestamp << std::put_time(std::gmtime(&time_t), "%Y-%m-%dT%H:%M:%S");
            bracket_log << "[" << timestamp.str() << "] " << msg << std::endl;
            bracket_log.close();
          }
        };
        
        update_jobs_with_result(self->state(), gg, theta, xs, uq, L, logger);
        
        auto bracket_end = std::chrono::high_resolution_clock::now();
        auto bracket_duration = std::chrono::duration_cast<std::chrono::microseconds>(bracket_end - bracket_start);
        double bracket_update_time_ms = bracket_duration.count() / 1000.0;
        
        self->println("Bracket updates processed in {:.3f}ms after result gg={}, theta={}, xs={:.2f}. Updated/Default jobs: {}/{}", 
                      bracket_update_time_ms, gg, theta, xs, self->state().updated_jobs.size(), self->state().default_jobs.size());
      }

      if (self->state().use_bracket) {
        auto strategy_it = self->state().strategy_map.find({gg, theta, xs});
        int strategy_used = (strategy_it != self->state().strategy_map.end()) ? strategy_it->second : 0;
        
        std::ofstream strategy_log("strategy_effectiveness_log.txt", std::ios::app);
        if (strategy_log.is_open()) {
          auto now = std::chrono::system_clock::now();
          auto time_t = std::chrono::system_clock::to_time_t(now);
          std::stringstream timestamp;
          timestamp << std::put_time(std::gmtime(&time_t), "%Y-%m-%dT%H:%M:%S");
          
          strategy_log << "[" << timestamp.str() << "] "
                      << "gg=" << gg << " "
                      << "theta=" << theta << " "
                      << "xs=" << xs << " "
                      << "strategy=" << strategy_used << " "
                      << "result=" << uq << " "
                      << "execution_time=" << job_runtime_seconds << "s "
                      << "quenched=" << (uq < 0.0 ? "true" : "false");
          
          if (job_found) {
            double usmin_calc = -1000.0 / std::pow(2.0, completed_job.n);
            strategy_log << " usmin=" << usmin_calc << " usmax=0.0";
          }
          
          strategy_log << std::endl;
          strategy_log.close();
        }
      }

      if (!self->state().updated_jobs.empty() || !self->state().default_jobs.empty()) {
        job_t job;
        bool has_job = false;
        bool is_updated_job = false;
        
        // First priority: get a job with updated brackets
        if (!self->state().updated_jobs.empty()) {
          job = self->state().updated_jobs.back();
          self->state().updated_jobs.pop_back();
          has_job = true;
          is_updated_job = true;
        } 
        // Second priority: get a job with default brackets
        else if (!self->state().default_jobs.empty()) {
          job = self->state().default_jobs.back();
          self->state().default_jobs.pop_back();
          has_job = true;
          is_updated_job = false;
        }
        
        if (has_job) {
          running_job_info job_info;
          job_info.job = job;
          job_info.worker = worker;
          job_info.start_time = std::chrono::high_resolution_clock::now();
          
          if (is_updated_job) {
            // Get the strategy that was applied to this job
            auto strategy_it = self->state().strategy_map.find({job.gg, job.theta, job.xs});
            job_info.current_strategy = (strategy_it != self->state().strategy_map.end()) ? strategy_it->second : 0;
            double usmin_calc = -1000.0 / std::pow(2.0, job.n);
            self->println("Reassigning updated job with strategy {} and bracket [{}, {}] for gg={}, theta={}, xs={}", 
                          job_info.current_strategy, usmin_calc, 0.0, job.gg, job.theta, job.xs);
          } else {
            job_info.current_strategy = 0; // Default strategy (no bracket update)
            double usmin_calc = -1000.0 / std::pow(2.0, job.n);
            self->println("Reassigning default job with bracket [{}, {}] for gg={}, theta={}, xs={}", 
                          usmin_calc, 0.0, job.gg, job.theta, job.xs);
          }
          
          self->state().running_jobs[worker] = job_info;
          ++self->state().active_jobs;
          self->anon_send(worker, std::move(job));
          self->println("Remaining jobs: {} (updated: {}, default: {})", 
                        self->state().updated_jobs.size() + self->state().default_jobs.size(),
                        self->state().updated_jobs.size(), self->state().default_jobs.size());
        }
      } else if (self->state().active_jobs > 0) {
        self->println("No more jobs available, waiting for active jobs to finish...");
        self->anon_send(worker, "quit");
      } else if (self->state().active_jobs == 0) {
        self->println("Total jobs done: {}", self->state().jobs_done);
        self->anon_send(worker, "quit");
        self->quit();
      }
    },
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
