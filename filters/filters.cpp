#include <cmath>
#include <thread>
#include <vector>

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"


namespace Filter {
    namespace Gauss { 
        //Dont touch this function, since it will impact the results for the comparison, even if it approximate to a similar value, it will not be exactly the same.
        void get_weights(int n, double* weights_out)
        {
            for (auto i { 0 }; i <= n; i++) { //change for itteration two.
                double x { static_cast<double>(i) * max_x / n};
                weights_out[i] = exp(-x * x * pi);
            }
        }
    }

    struct Coefficients_Blur{
        int segment_starts{0};
        int segment_stops{0};
        float xAxis_size{0};
        float yAxis_size{0};
        unsigned radius{0};
        Matrix* m;
        Matrix* writeMatrix;
    };

    void* getBlurCoefficients_threads(void* arg) {
    Coefficients_Blur* passer = static_cast<Coefficients_Blur*>(arg);
    Matrix scratch{PPM::max_dimension};
    double w[Gauss::max_radius]{};
    Gauss::get_weights(passer->radius, w);

    for (auto y{0}; y < passer->yAxis_size; y++) {
        for (auto x{passer->segment_starts}; x < passer->segment_stops; x++) {
            auto r{w[0] * passer->m->r(x, y)}, g{w[0] * passer->m->g(x, y)}, b{w[0] * passer->m->b(x, y)}, n{w[0]};
            for (auto wi{1}; wi <= passer->radius; wi++) {
                auto wc{w[wi]};
                auto x2{x - wi};
                if (x2 >= 0) {
                    b += wc * passer->m->b(x2, y);
                    g += wc * passer->m->g(x2, y);
                    r += wc * passer->m->r(x2, y);
                    n += wc;
                }
                x2 = x + wi;
                if (x2 < passer->xAxis_size) {
                    b += wc * passer->m->b(x2, y);
                    g += wc * passer->m->g(x2, y);
                    r += wc * passer->m->r(x2, y);
                    n += wc;
                }
            }
            scratch.b(x, y) = b / n;
            scratch.g(x, y) = g / n;
            scratch.r(x, y) = r / n;
        }
    }

    for (auto y{0}; y < passer->yAxis_size; y++) {
        for (auto x{passer->segment_starts}; x < passer->segment_stops; x++) {
            auto r{w[0] * scratch.r(x, y)}, g{w[0] * scratch.g(x, y)}, b{w[0] * scratch.b(x, y)}, n{w[0]};
            for (auto wi{1}; wi <= passer->radius; wi++) {
                auto y2{y - wi};
                auto wc{w[wi]};
                if (y2 >= 0) {
                    b += wc * scratch.b(x, y2);
                    g += wc * scratch.g(x, y2);
                    r += wc * scratch.r(x, y2);
                    n += wc;
                }
                y2 = y + wi;
                if (y2 < passer->yAxis_size) {
                    b += wc * scratch.b(x, y2);
                    g += wc * scratch.g(x, y2);
                    r += wc * scratch.r(x, y2);
                    n += wc;
                }
            }
            passer->writeMatrix->b(x, y) = b / n;
            passer->writeMatrix->g(x, y) = g / n;
            passer->writeMatrix->r(x, y) = r / n;
        }
    }

    return nullptr;
}

    void updateIntervals(const float xAxis_size, const int number_of_threads, Coefficients_Blur* intervals);
    void blur_adjustIntervals(const float xAxis_size, const unsigned newDiff, const int number_of_threads, Coefficients_Blur* intervals);
    void createWorkerThreads(const int number_of_threads, pthread_t* workerThreads, Coefficients_Blur* intervals);

    Matrix blur(Matrix m, const int radius,const int number_of_threads)
    {
        Matrix dst { m };
        Matrix writeMatrix { m };
        const float xAxis_size = dst.get_x_size();
        const float yAxis_size = dst.get_y_size();
        if(number_of_threads==1){
            getBlurWithSingleThread(m, radius, dst, writeMatrix);;
        }
        else{
            Coefficients_Blur* intervals = new Coefficients_Blur[number_of_threads];
            int offsetX = xAxis_size / number_of_threads;
            auto num_Iterations = number_of_threads;
            auto* pointed_intervals = intervals;
            auto ind_position = 0;
            do {
                pointed_intervals->xAxis_size = xAxis_size;
                pointed_intervals->yAxis_size = yAxis_size;
                pointed_intervals->segment_starts = offsetX * ind_position;
                pointed_intervals->segment_stops = offsetX * ind_position + offsetX;
                pointed_intervals->radius = radius;
                pointed_intervals->m = &dst;
                pointed_intervals->writeMatrix = &writeMatrix;
                pointed_intervals++;
                ind_position++;

            } while (--num_Iterations > 0);

            unsigned newDiff = (xAxis_size > intervals[number_of_threads - 1].segment_stops) ? (xAxis_size - intervals[number_of_threads - 1].segment_stops) : 0;

            blur_adjustIntervals(xAxis_size, newDiff, number_of_threads, intervals);
            pthread_t workerThreads[number_of_threads];
            createWorkerThreads(number_of_threads, workerThreads, intervals);
            for (auto i = 0; i < number_of_threads; i +=2){
                if (pthread_join(workerThreads[i], nullptr) != 0) {
                        // tasks to do on error, if any
                }
                if (i + 1 < number_of_threads) {
                    if (pthread_join(workerThreads[i + 1], nullptr) != 0) {
                        // tasks to do on error, if any
                    }
                }
            }
        }
        return writeMatrix;
    }
    void getBlurWithSingleThread(const Matrix& m, const int radius, Matrix& dst, Matrix& writeMatrix)
    {
    Coefficients_Blur passer;
    passer.segment_starts = 0;
    passer.segment_stops = m.get_x_size();
    passer.xAxis_size = m.get_x_size();
    passer.yAxis_size = m.get_y_size();
    passer.radius = radius;
    passer.m = const_cast<Matrix*>(&m);
    passer.writeMatrix = &writeMatrix;

    getBlurCoefficients_threads(static_cast<void*>(&passer));
    }
    void updateIntervals(const float xAxis_size, const int number_of_threads, Coefficients_Blur* intervals)
    {
        const int offsetX = xAxis_size / number_of_threads;

        for (int i = 0; i < number_of_threads; i++) {
            intervals[i].xAxis_size = xAxis_size;
            intervals[i].segment_starts = offsetX * i;
            intervals[i].segment_stops = offsetX * i + offsetX;
        }
    }
    void blur_adjustIntervals(const float xAxis_size, const unsigned newDiff, const int number_of_threads, Coefficients_Blur* intervals)
    {
        auto i = 0;
        auto* intervalsPtr = intervals;
        auto numIterations = newDiff;
        while (numIterations--) {
            intervalsPtr->segment_starts += 1;
            intervalsPtr->segment_stops += 1;
            intervalsPtr++;
        }

        intervals[0].segment_starts = (newDiff > 0) ? 0 : intervals[0].segment_starts;
        intervals[number_of_threads - 1].segment_stops = (newDiff > 0) ? xAxis_size : intervals[number_of_threads - 1].segment_stops;
    }
    void createWorkerThreads(const int number_of_threads, pthread_t* workerNodeParallelThread, Coefficients_Blur* intervals)
    {
        for (auto i = 0; i < number_of_threads; i += 2) {
            if (pthread_create(&workerNodeParallelThread[i], nullptr, getBlurCoefficients_threads, &intervals[i]) != 0) {
                // tasks to do on error, if any
            }
            if (i + 1 < number_of_threads) {
                if (pthread_create(&workerNodeParallelThread[i + 1], nullptr, getBlurCoefficients_threads, &intervals[i + 1]) != 0) {
                    // tasks to do on error, if any
                }
            }
        }
    }


    struct Threshold_Passer{
        int pixel_first{0};
        int pixel_last{0};
        int pixel_sum{0};
        Matrix* m;
    };

    void processData(unsigned start, unsigned stop, unsigned sum, Matrix* m);

    void* thresholdThreadLimit(void* arg) {
        Threshold_Passer* data_pipe_line = static_cast<Threshold_Passer*>(arg);
        unsigned start = data_pipe_line->pixel_first;
        unsigned stop = data_pipe_line->pixel_last;
        unsigned sum = data_pipe_line->pixel_sum;
        Matrix* m = data_pipe_line->m;
        processData(start, stop, sum, m);
        pthread_exit(NULL);
    }

    void* generateThreshold(unsigned start, unsigned stop, unsigned sum, Matrix* m) {
        processData(start, stop, sum, m);
        return 0;
    }
    void processData(unsigned start, unsigned stop, unsigned sum, Matrix* m) {
        unsigned psum {};
        unsigned int i = start;
        while (i < stop) {
            psum = m->r(i, 0) + m->g(i, 0) + m->b(i, 0);
            int pixelValue = (sum > psum) ? 0 : 255;
            m->r(i, 0) = m->g(i, 0) = m->b(i, 0) = pixelValue;
            i++;
        }
    }

    void adjustIntervalsInThreads(const float nump, const unsigned newDiff, const int number_of_threads, Threshold_Passer* intervals);
    void createParallelThreads(const int number_of_threads, pthread_t* workerNodeParallelThread, Threshold_Passer* intervals);

    Matrix threshold(Matrix m, const int number_of_threads)
    {
        auto dst = m;
        unsigned sum = 0;
        unsigned nump = dst.get_x_size() * dst.get_y_size();
        getPixelsSum(m, sum, nump);
        if(number_of_threads==1){
            sum /= nump;
            generateThreshold(0,nump,sum,&dst);
        }
        else{
            //Compute ranges.
            Threshold_Passer* intervals = new Threshold_Passer[number_of_threads];
            int offset = nump/number_of_threads;
            auto num_Iterations = number_of_threads;
            auto* pointed_intervals = intervals;
            auto ind_position = 0;
            do {
                pointed_intervals->pixel_first = offset * ind_position;
                pointed_intervals->pixel_last = offset * ind_position + offset;
                pointed_intervals->pixel_sum = sum/nump;
                pointed_intervals->m = &dst;
                pointed_intervals++;
                ind_position++;

            } while (--num_Iterations > 0);
            unsigned newDiff = (nump > intervals[number_of_threads - 1].pixel_last) ? (nump - intervals[number_of_threads - 1].pixel_last) : 0;

            adjustIntervalsInThreads(nump, newDiff, number_of_threads, intervals);
            pthread_t workerNodeParallelThread[number_of_threads];
            createParallelThreads(number_of_threads, workerNodeParallelThread, intervals);
            for (int i = 0; i < number_of_threads; i++) {
                if (pthread_join(workerNodeParallelThread[i], nullptr) != 0) {
                        // tasks to do on error, if any
                    }
                    if (i + 1 < number_of_threads) {
                        if (pthread_join(workerNodeParallelThread[i + 1], nullptr) != 0) {
                            // tasks to do on error, if any
                        }
                    }
            }
        }
        return dst;
    }
    void getPixelsSum(const Matrix& m, unsigned& sum, unsigned nump) {
        for (auto i = 0; i < nump; i++) {
            auto b = m.b(i, 0);
            auto g = m.g(i, 0);
            auto r = m.r(i, 0);
            sum += r + g + b;
        }
    }
    void adjustIntervalsInThreads(const float nump, const unsigned newDiff, const int number_of_threads, Threshold_Passer* intervals)
    {
        auto i = 0;
        auto numIterations = newDiff;
        auto* intervalsPtr = intervals;
        while (numIterations--) {
            intervalsPtr->pixel_first += 1;
            intervalsPtr->pixel_last += 1;
            intervalsPtr++;
        }

        intervals[0].pixel_first = (newDiff > 0) ? 0 : intervals[0].pixel_first;
        intervals[number_of_threads - 1].pixel_last = (newDiff > 0) ? nump : intervals[number_of_threads - 1].pixel_last;
    }
    void createParallelThreads(const int number_of_threads, pthread_t* workerNodeParallelThread, Threshold_Passer* intervals)
    {
        for (auto i = 0; i < number_of_threads; i += 2) {
            if (pthread_create(&workerNodeParallelThread[i], nullptr, thresholdThreadLimit, &intervals[i]) != 0) {
            // tasks to do on error, if any
            }
            if (i + 1 < number_of_threads) {
                if (pthread_create(&workerNodeParallelThread[i + 1], nullptr, thresholdThreadLimit, &intervals[i + 1]) != 0) {
                    // tasks to do on error, if any
                }
            }
        }
    }
}


