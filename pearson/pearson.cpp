#include <bits/stdc++.h>
#include <cstdlib>
#include <iostream>
#include <pthread.h>

#include "analysis.hpp"
#include "dataset.hpp"

using namespace Dataset;

void devide_tasks_among_threads_and_store_result(int argc, char const*argv[]);

int main(int argc, char const* argv[])
{
    if(argc>=4){
        devide_tasks_among_threads_and_store_result(argc,argv);
        return 0;
    }
    else if (argc < 4) {
        char temp_warning[1024];
        snprintf(temp_warning,1024,"Usage: %s <dataset> <outfile> <max_number_of_threads>\n",argv[0]);
        perror((const char*)temp_warning);
        return 1;
    }
}

void devide_tasks_among_threads_and_store_result(int argc, char const*argv[]){
    int max_allowed_threads = atoi(argv[3]);

    std::vector<double> final_results;
    pthread_t par_threads[max_allowed_threads];
    struct thread_states threads_state[max_allowed_threads];
    
    auto datasets { read(argv[1]) };
    unsigned int set_size = datasets[0].get_size();
    // spliting the tasks into threads and do correlation_coefficients
    for(int i = 0; i < max_allowed_threads; i++){
        threads_state[i].serial_number = i;
        threads_state[i].splitted_task_size = set_size;
        threads_state[i].max_allowed_threads = max_allowed_threads;
        threads_state[i].coefficients = datasets;
        pthread_create(&par_threads[i], NULL, Analysis::correlation_coefficients, (void *)&threads_state[i]);
    }
    void *status;
    for(int i = 0; i < max_allowed_threads; i++)
        pthread_join(par_threads[i], &status); // Let's wait for all the threads to finish their tasks

    // merging each thread's into one
    for(int i = 0; i < max_allowed_threads; i++)
        final_results.insert(final_results.end(), threads_state[i].par_splitted_results.begin(), threads_state[i].par_splitted_results.end());

    write(final_results, argv[2]);
    pthread_exit(NULL);
}