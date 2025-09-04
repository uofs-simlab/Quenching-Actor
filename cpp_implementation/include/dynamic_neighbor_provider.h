#ifndef DYNAMIC_NEIGHBOR_PROVIDER_H
#define DYNAMIC_NEIGHBOR_PROVIDER_H

#include "bracket_optimizer.h"
#include <unordered_map>
#include <unordered_set>

namespace bracket_optimizer {

    struct job_node {
        job_t job;
        enum status_t { pending, running, done } status;
        double uq; 
        std::unordered_set<std::tuple<int, float, float>> neighbors;
        
        job_node() = default;
        
        job_node(const job_t& j, status_t s, double u, const std::unordered_set<std::tuple<int, float, float>>& n = {}) 
            : job(j), status(s), uq(u), neighbors(n) {}
    };

    class DynamicNeighborProvider : public NeighborProvider {
    private:
        const std::unordered_map<std::tuple<int, float, float>, job_node>* job_graph_;
        
    public:
        explicit DynamicNeighborProvider(
            const std::unordered_map<std::tuple<int, float, float>, job_node>* job_graph)
            : job_graph_(job_graph) {}

        std::vector<std::pair<std::tuple<int, float, float>, double>> 
        get_completed_neighbors(const std::tuple<int, float, float>& job_key) override {
            
            std::vector<std::pair<std::tuple<int, float, float>, double>> completed_neighbors;
            
            if (!job_graph_) return completed_neighbors;
            
            auto job_it = job_graph_->find(job_key);
            if (job_it == job_graph_->end()) {
                return completed_neighbors;
            }
            
            const auto& node = job_it->second;
            
            for (const auto& neighbor_key : node.neighbors) {
                auto neighbor_it = job_graph_->find(neighbor_key);
                if (neighbor_it != job_graph_->end() && 
                    neighbor_it->second.status == job_node::done) {
                    completed_neighbors.emplace_back(neighbor_key, neighbor_it->second.uq);
                }
            }
            
            return completed_neighbors;
        }
    };

} // namespace bracket_optimizer

#endif // DYNAMIC_NEIGHBOR_PROVIDER_H
