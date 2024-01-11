#include <vector>
#include "vector.hpp"


#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

    struct thread_states {
        int serial_number;
        int splitted_task_size;
        int max_allowed_threads;
        std::vector<double> par_splitted_results;
        std::vector<Vector> coefficients;
    };

    namespace Analysis {
        void itterate_over_segmented_matrix(struct thread_states *local_thread_state,int starting_position, int ending_position);
        void* correlation_coefficients(void *datasets);
        double pearson(Vector vec1, Vector vec2);
    };

#endif
