#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <list>

#include "analysis.hpp"

namespace Analysis {

    void* correlation_coefficients(void *datasets)
    {
        struct thread_states *local_thread_state = (struct thread_states*) datasets;
        std::vector<double> coefficients {};
        int max_number_of_threads = local_thread_state->max_allowed_threads;
        int splited_payload_width = local_thread_state->splitted_task_size / max_number_of_threads;
        int starting_position = local_thread_state->serial_number * splited_payload_width;
        int ending_position = starting_position + splited_payload_width;
        itterate_over_segmented_matrix(local_thread_state,starting_position,ending_position);
        pthread_exit((void*)local_thread_state); // signal the thread execution completed to the caller
    }

    double pearson(Vector vec1, Vector vec2)
    {
        auto x_mean { vec1.mean() };
        auto x_mm { vec1 - x_mean };
        auto x_mag { x_mm.magnitude() };
        auto x_mm_over_x_mag { x_mm / x_mag };

        auto y_mean { vec2.mean() };
        auto y_mm { vec2 - y_mean };
        auto y_mag { y_mm.magnitude() };
        auto y_mm_over_y_mag { y_mm / y_mag };

        auto r { x_mm_over_x_mag.dot(y_mm_over_y_mag) };

        return std::max(std::min(r, 1.0), -1.0);
    }

    void itterate_over_segmented_matrix(struct thread_states *local_thread_state,int starting_position, int ending_position){
        for (int sample1 { starting_position }; sample1 < ending_position; sample1++) {
            for (int sample2 { sample1 + 1 }; sample2 < local_thread_state->splitted_task_size; sample2++) {
                double corr { pearson(local_thread_state->coefficients[sample1], local_thread_state->coefficients[sample2]) };
                local_thread_state->par_splitted_results.push_back(corr);
            }
        }
    }

};
