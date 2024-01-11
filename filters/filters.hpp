/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"

#if !defined(FILTERS_HPP)
#define FILTERS_HPP

namespace Filter {

    namespace Gauss {
        constexpr unsigned max_radius { 1000 };
        constexpr float max_x { 1.33 };
        constexpr float pi { 3.14159 };
        
        void get_weights(int n, double* weights_out);
        void* thresholdThreadLimit(struct test);
        void* generateThreshold(unsigned start, unsigned stop,unsigned sum, Matrix m);
        void* getBlurCoefficients(unsigned segment_stops,unsigned stopY,unsigned radius, Matrix* dst);
        void* getBlurCoefficientsThread(void* arg);
    }
    void getBlurWithSingleThread(const Matrix& m, const int radius, Matrix& dst, Matrix& writeMatrix);
    void getBlurWithMultiThread(const Matrix& m, const int radius, const int number_of_threads, Matrix& dst, Matrix& writeMatrix);
    Matrix blur(Matrix m, const int radius, const int number_of_threads);
    void getPixelsSum(const Matrix& m, unsigned& sum, unsigned nump);
    Matrix threshold(Matrix m, const int number_of_threads);
    struct DataPasser;
}

#endif
