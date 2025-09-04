#ifndef STATIC_NEIGHBOR_PROVIDER_H
#define STATIC_NEIGHBOR_PROVIDER_H

#include "bracket_optimizer.h"

namespace bracket_optimizer {

    class StaticNeighborProvider : public NeighborProvider {
    public:
        std::vector<std::pair<std::tuple<int, float, float>, double>> 
        get_completed_neighbors(const std::tuple<int, float, float>& job_key) override {
            // Static system doesn't track neighbors during execution
            return {};
        }
    };

} // namespace bracket_optimizer

#endif // STATIC_NEIGHBOR_PROVIDER_H
