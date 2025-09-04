#ifndef WAVE_LOADER_H
#define WAVE_LOADER_H

#include "fitzhugh_nagumo.h"
#include <string>

/**
 * @brief Utility functions for loading wave data from files
 */
class WaveLoader {
public:
    static Matrix load_wave_data(const std::string& filename);

    static Vector load_parameters(const std::string& filename);

    static std::tuple<Matrix, Vector, double, Vector> load_and_center_simple(
        const std::string& Ufile, 
        const std::string& Pfile, 
        const FitzHughNagumoSolver& solver);

    static States create_states_from_files(
        const std::string& fastU, const std::string& fastP,
        const std::string& slowU, const std::string& slowP,
        const FitzHughNagumoSolver& solver);
};

#endif // WAVE_LOADER_H
