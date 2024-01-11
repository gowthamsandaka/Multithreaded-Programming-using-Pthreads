#include <string>
#include <vector>
#include "vector.hpp"

#ifndef DATASET_HPP
#define DATASET_HPP

    namespace Dataset {
        std::vector<Vector> read(std::string filename);
        void write(std::vector<double> data, std::string filename);
    };

#endif
