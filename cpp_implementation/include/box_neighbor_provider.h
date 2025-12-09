#ifndef BOX_NEIGHBOR_PROVIDER_H
#define BOX_NEIGHBOR_PROVIDER_H

#include "bracket_optimizer.h"
#include <unordered_map>
#include <unordered_set>
#include <tuple>

namespace bracket_optimizer {

    /**
     * @brief Neighbor provider for box-based dynamic refinement systems
     * 
     * Uses spatial proximity within the same gg level to find neighbors.
     * Points within the same box or nearby boxes are considered neighbors.
     */
    class BoxNeighborProvider : public NeighborProvider {
    private:
        // Server state references
        const std::unordered_map<std::tuple<int, float, float>, std::pair<int, double>>* job_results_;
        const std::unordered_map<std::tuple<int, float, float>, std::unordered_set<std::string>>* point_to_box_map_;
        const std::unordered_map<std::string, void*>* active_boxes_; // void* to avoid circular dependency
        
        // Configuration
        float spatial_tolerance_theta_;
        float spatial_tolerance_xs_;
        
    public:
        BoxNeighborProvider(
            const std::unordered_map<std::tuple<int, float, float>, std::pair<int, double>>* job_results,
            const std::unordered_map<std::tuple<int, float, float>, std::unordered_set<std::string>>* point_to_box_map,
            const std::unordered_map<std::string, void*>* active_boxes,
            float theta_tolerance = 5.0f,  // Spatial search radius for theta
            float xs_tolerance = 5.0f       // Spatial search radius for xs
        ) : job_results_(job_results), 
            point_to_box_map_(point_to_box_map),
            active_boxes_(active_boxes),
            spatial_tolerance_theta_(theta_tolerance),
            spatial_tolerance_xs_(xs_tolerance) {}
        
        std::vector<std::pair<std::tuple<int, float, float>, double>> 
        get_completed_neighbors(const std::tuple<int, float, float>& job_key) override {
            
            std::vector<std::pair<std::tuple<int, float, float>, double>> neighbors;
            
            auto [target_gg, target_theta, target_xs] = job_key;
            
            // Strategy 1: Find neighbors within the same boxes
            auto box_it = point_to_box_map_->find(job_key);
            if (box_it != point_to_box_map_->end()) {
                for (const std::string& box_id : box_it->second) {
                    // Find all other points in the same box(es)
                    for (const auto& [point_key, box_set] : *point_to_box_map_) {
                        if (box_set.count(box_id) > 0 && point_key != job_key) {
                            auto [point_gg, point_theta, point_xs] = point_key;
                            
                            // Only consider same gg level
                            if (point_gg == target_gg) {
                                auto result_it = job_results_->find(point_key);
                                if (result_it != job_results_->end() && result_it->second.first == 2) { // done = 2
                                    neighbors.emplace_back(point_key, result_it->second.second);
                                }
                            }
                        }
                    }
                }
            }
            
            // Strategy 2: Find spatially close neighbors (within tolerance)
            for (const auto& [point_key, result_pair] : *job_results_) {
                auto [point_gg, point_theta, point_xs] = point_key;
                
                // Only consider same gg level and completed jobs
                if (point_gg == target_gg && result_pair.first == 2 && point_key != job_key) { // done = 2
                    float theta_dist = std::abs(point_theta - target_theta);
                    float xs_dist = std::abs(point_xs - target_xs);
                    
                    // Check if within spatial tolerance
                    if (theta_dist <= spatial_tolerance_theta_ && xs_dist <= spatial_tolerance_xs_) {
                        // Avoid duplicates from Strategy 1
                        bool already_added = false;
                        for (const auto& [existing_key, existing_uq] : neighbors) {
                            if (existing_key == point_key) {
                                already_added = true;
                                break;
                            }
                        }
                        if (!already_added) {
                            neighbors.emplace_back(point_key, result_pair.second);
                        }
                    }
                }
            }
            
            return neighbors;
        }
    };

} // namespace bracket_optimizer

#endif // BOX_NEIGHBOR_PROVIDER_H