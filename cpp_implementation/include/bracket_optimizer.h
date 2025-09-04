#ifndef BRACKET_OPTIMIZER_H
#define BRACKET_OPTIMIZER_H

#include <unordered_map>
#include <vector>
#include <tuple>
#include <memory>
#include "tuple_hash.h"  
#include "job_structures.h"

namespace bracket_optimizer {

    using uq_cache_t = std::unordered_map<std::tuple<int, float>, std::vector<result_entry>>;

    class NeighborProvider {
    public:
        virtual ~NeighborProvider() = default;
        virtual std::vector<std::pair<std::tuple<int, float, float>, double>> 
            get_completed_neighbors(const std::tuple<int, float, float>& job_key) = 0;
    };


    struct BracketConfig {
        bool enable_neighbor_strategy = true;   // Use neighbor-based optimization
        bool enable_cache_strategy = true;      // Use cache-based optimization
        double trust_margin = 1.5;              // Trust margin multiplier for neighbor results
        double min_bracket_ratio = 1000.0;      // Minimum ratio for bracket calculation
        int max_gg_search = -1;                 // Maximum gg value to search in cache (-1 = auto-detect from available cache levels)
        
        BracketConfig() = default;
    };

    /**
     * @brief Initialize job bracket parameters based on cache and neighbor information
     * 
     * This function optimizes the initial bracket [usmin, usmax] for the bisection solver
     * by leveraging previous results from neighboring jobs and cached results.
     * 
     * @param job Reference to the job to optimize (n parameter will be modified)
     * @param uq_cache Cache of previous uq results indexed by (gg, theta)
     * @param neighbor_provider Optional provider for getting neighbor information
     * @param config Configuration parameters for optimization
     * @return true if bracket was successfully optimized, false if using default
     */
    bool initialize_bracket_from_cache(
        job_t& job,
        const uq_cache_t& uq_cache,
        std::shared_ptr<NeighborProvider> neighbor_provider = nullptr,
        const BracketConfig& config = BracketConfig()
    );

    /**
     * @brief Calculate the n parameter from a uq value
     * @param uq The uq value to convert
     * @param config Configuration parameters
     * @return The calculated n parameter
     */
    int calculate_n_from_uq(double uq, const BracketConfig& config = BracketConfig());

    /**
     * @brief Convert n parameter to usmin value
     * @param n The n parameter
     * @param config Configuration parameters
     * @return The corresponding usmin value
     */
    double n_to_usmin(int n, const BracketConfig& config = BracketConfig());

} // namespace bracket_optimizer

#endif // BRACKET_OPTIMIZER_H
