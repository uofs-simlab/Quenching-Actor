#include "bracket_optimizer.h"
#include "job_structures.h"
#include <cmath>
#include <algorithm>

namespace bracket_optimizer {

    int calculate_n_from_uq(double uq, const BracketConfig& config) {
        if (uq >= 0.0) {
            return 0; // Default for non-quenching results: usmin = -1000
        }
        
        double adjusted_uq = uq * config.trust_margin;
        if (adjusted_uq < 0.0) {
            double abs_adjusted_uq = std::abs(adjusted_uq);
            if (abs_adjusted_uq > 0.0) {
                return static_cast<int>(std::floor(std::log2(config.min_bracket_ratio / abs_adjusted_uq)));
            }
        }
        return 0; // Fallback: usmin = -1000
    }

    double n_to_usmin(int n, const BracketConfig& config) {
        return -config.min_bracket_ratio / std::pow(2.0, n);
    }

    bool initialize_bracket_from_cache(
        job_t& job,
        const uq_cache_t& uq_cache,
        std::shared_ptr<NeighborProvider> neighbor_provider,
        const BracketConfig& config
    ) {
        bool bracket_updated = false;

        // STRATEGY 1: Cross-GG optimization with same theta and exact xs match
        if (!bracket_updated && config.enable_cache_strategy) {
            int search_limit = config.max_gg_search;
            if (search_limit == -1) {
                search_limit = job.gg;
                for (const auto& [cache_key, cache_entries] : uq_cache) {
                    int cache_gg = std::get<0>(cache_key);
                    if (cache_gg > search_limit) {
                        search_limit = cache_gg;
                    }
                }
            }
            
            for (int g2 = job.gg + 1; g2 <= search_limit; ++g2) {
                auto cache_key = std::make_tuple(g2, job.theta);
                auto it = uq_cache.find(cache_key);
                if (it != uq_cache.end()) {
                    for (const auto& entry : it->second) {
                        if (entry.uq < 0.0 && std::abs(entry.xs - job.xs) < 1e-4f) {
                            job.n = calculate_n_from_uq(entry.uq, config);
                            bracket_updated = true;
                            break; // Use first match
                        }
                    }
                    if (bracket_updated) break;
                }
            }
        }

        // STRATEGY 2: Spatial neighbors within same gg level
        if (!bracket_updated && config.enable_neighbor_strategy && neighbor_provider) {
            try {
                auto job_key = std::make_tuple(job.gg, job.theta, job.xs);
                auto completed_neighbors = neighbor_provider->get_completed_neighbors(job_key);
                
                for (const auto& [neighbor_key, neighbor_uq] : completed_neighbors) {
                    if (neighbor_uq < 0.0) {
                        job.n = calculate_n_from_uq(neighbor_uq, config);
                        bracket_updated = true;
                        break; // Use first neighbor match
                    }
                }
            } catch (...) {
                // Silently continue if neighbor lookup fails
            }
        }

        return bracket_updated;
    }

} // namespace bracket_optimizer
