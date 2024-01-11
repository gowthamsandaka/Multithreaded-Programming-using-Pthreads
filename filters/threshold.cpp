/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char const* argv[])
{
    unsigned int number_of_threads = (argc == 4 && argv[3]) ? static_cast<unsigned>(std::stoul(argv[3])) : 0;
    if (argc != 4 || number_of_threads<1) {
        std::cerr << "Usage: " << argv[0] << " [infile] [outfile] [Threads 1-N]" << std::endl;
        std::exit(1);
    }

    PPM::Reader reader {};
    PPM::Writer writer {};

    auto m { reader(argv[1]) };
    auto thresholded { Filter::threshold(m,number_of_threads) };

    writer(thresholded, argv[2]);

    return 0;
}
