/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"
#include "ppm.hpp"
#include "filters.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char const* argv[])
{
    unsigned int number_of_threads = (argc == 5 && argv[4]) ? static_cast<unsigned>(std::stoul(argv[4])) : 0;
    if (argc != 5 || number_of_threads<1) {
        std::cerr << "Usage: " << argv[0] << " [radius] [infile] [outfile] [Threads 1-n]" << std::endl;
        std::exit(1);
    }

    PPM::Reader reader {};
    PPM::Writer writer {};

    auto m { reader(argv[2]) };// Läser in en bild
    auto radius { static_cast<unsigned>(std::stoul(argv[1])) };//sparar radien.

    auto blurred { Filter::blur(m, radius,number_of_threads) };//applicerar filtret
    writer(blurred, argv[3]);//Printar filtret

    return 0;
}
