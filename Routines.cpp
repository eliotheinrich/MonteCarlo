#include "MonteCarlo.cpp"

#define BOLTZMANN_CONSTANT 0.08617

template <class ModelType>
void exchange(std::vector<MonteCarlo<ModelType>*> *models, std::vector<float> *T) {
    MonteCarlo<ModelType> *model_buffer;
    float r;
    float dE; float dB;
    for (int i = 0; i < models->size()-1; i++) {
        r = float((*models)[i]->r())/float(RAND_MAX);
        dE = (*models)[i]->energy - (*models)[i+1]->energy;
        dB = 1./(*T)[i] - 1./(*T)[i+1];
        if (r < exp(dE*dB)) {
            model_buffer = (*models)[i];
            (*models)[i] = (*models)[i+1];
            (*models)[i+1] = model_buffer;
        }
    }
}

template <class ModelType, class dtype>
std::vector<std::vector<std::vector<dtype>>> sample_r(std::vector<dtype> sampling_func(ModelType*), 
                                                      ModelType *model, 
                                                      std::vector<float> &T, 
                                                      unsigned long long steps_per_run, 
                                                      unsigned long long num_samples, 
                                                      unsigned long long steps_per_sample,
                                                      int num_threads) {

    int resolution = T.size();
    int dtype_size = sampling_func(model).size();

    std::vector<MonteCarlo<ModelType>*> models(resolution);

    for (int i = 0; i < resolution; i++) {
        models[i] = new MonteCarlo<ModelType>(model->clone());
        models[i]->model->randomize_spins();
    }

    std::vector<std::vector<std::vector<dtype>>> arr = std::vector<std::vector<std::vector<dtype>>>(resolution, 
                                                                   std::vector<std::vector<dtype>>(dtype_size,
                                                                               std::vector<dtype>(2)));

    ctpl::thread_pool threads(num_threads);

    std::vector<std::future<void>> results(resolution);

    auto do_steps = [steps_per_run](int id, MonteCarlo<ModelType> *m, float T) {
        run_MC(m, steps_per_run, "trig", 3*T, T);
        //m->steps(steps_per_run, T);
    };

    // Do initial steps
    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(do_steps, models[i], T[i]);
    }

    // Join threads
    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }

    const int buffer_size = 1000;
    auto take_samples = [num_samples, steps_per_sample, dtype_size, sampling_func](int id, int i,
                                                                       MonteCarlo<ModelType> *m, float T, 
                                                                       std::vector<std::vector<dtype>> *arr_i) {

        std::vector<std::vector<dtype>> buffer = std::vector<std::vector<dtype>>(buffer_size, 
                                                             std::vector<dtype>(dtype_size));
        int n = 0; int k = 0;
        for (int n = 1; n < num_samples+1; n++) {
            for (int j = 0; j < dtype_size; j++) {
                auto sample = sampling_func(m->model);
                (*arr_i)[j][0] = ((n-1)*(*arr_i)[j][0] + sample[j])/n;
                (*arr_i)[j][1] = ((n-1)*(*arr_i)[j][1] + pow(sample[j], 2))/n;
            }
            m->steps(steps_per_sample, T);
        }
    };

    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(take_samples, i, models[i], T[i], &arr[i]);
    }

    // Join threads
    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }

    std::cout << models[0]->nsteps << std::endl;

    return arr;
}

template <class ModelType, class A>
std::vector<std::vector<A>> sample_r_old(A sampling_func(ModelType*), ModelType *model, std::vector<float> *T, 
                                                        unsigned long long num_runs,
                                                        unsigned long long steps_per_run, 
                                                        unsigned long long num_samples, 
                                                        unsigned long long steps_per_sample,
                                                        int num_threads) {

    int resolution = T->size();

    std::vector<MonteCarlo<ModelType>*> models(resolution);

    for (int i = 0; i < resolution; i++) {
        models[i] = new MonteCarlo<ModelType>(model->clone());
        models[i]->model->randomize_spins();
    }

    std::vector<std::vector<A>> arr = std::vector<std::vector<A>>(resolution, std::vector<A>(num_runs*num_samples));

    ctpl::thread_pool threads(num_threads);

    std::vector<std::future<void>> results(resolution);

    auto do_steps = [steps_per_run](int id, MonteCarlo<ModelType> *m, float T) {
        m->steps(steps_per_run, T);
    };

    // Do initial steps
    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(do_steps, models[i], (*T)[i]);
    }

    // Join threads
    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }

    auto take_samples = [num_samples, steps_per_sample, num_runs, steps_per_run, sampling_func](int id, int i,
                                                                                    MonteCarlo<ModelType> *m, float T, 
                                                                                    std::vector<A> *arr_i) {
        for (int n = 0; n < num_runs; n++) {
            m->model->randomize_spins();
            m->steps(steps_per_run, T);
            for (int j = 0; j < num_samples; j++) {
                A sample = sampling_func(m->model);
                (*arr_i)[n*num_samples + j] = sample;
                m->steps(steps_per_sample, T);
            }
        }
    };


    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(take_samples, i, models[i], (*T)[i], &arr[i]);
    }

    // Join threads
    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }

    return arr;
}

template <class dtype>
std::vector<std::vector<std::vector<dtype>>> summary_statistics(std::vector<std::vector<std::vector<dtype>>> *data) {
    int resolution = (*data).size();
    int num_samples = (*data)[0].size();
    int dtype_size = (*data)[0][0].size();

    std::vector<std::vector<std::vector<dtype>>> arr(resolution, std::vector<std::vector<dtype>>(2, std::vector<dtype>(dtype_size)));

    for (int i = 0; i < resolution; i++) {
        for (int k = 0; k < dtype_size; k++) {
            for (int n = 0; n < num_samples; n++) {
                arr[i][0][k] += (*data)[i][n][k]/num_samples;
            }
        }
    }

    for (int i = 0; i < resolution; i++) {
        for (int k = 0; k < dtype_size; k++) {
            for (int n = 0; n < num_samples; n++) {
                arr[i][1][k] += pow((*data)[i][n][k] - arr[i][0][k], 2)/num_samples;
            }
            arr[i][1][k] = sqrt(arr[i][1][k]);
        }
    }

    return arr;
}

template <class dtype>
void write_data(std::vector<std::vector<std::vector<dtype>>> *data, std::vector<float> *T, std::string filename, std::string header = "") {
    int resolution = (*data).size();
    int dtype_size = (*data)[0].size();

    std::ofstream output_file(filename);

    // Write header
    output_file << resolution << "\t" << header << "\n";

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t";
        for (int k = 0; k < dtype_size; k++) {
            output_file << (*data)[i][k][0] << ", " << (*data)[i][k][1];
            if (k != dtype_size - 1) { output_file << "\t"; }
        }
        output_file << "\n";
    }

    output_file.close();
}

template <class dtype>
void write_samples(std::vector<std::vector<std::vector<dtype>>> *data, std::vector<float> *T, std::string filename, std::string header = "") {
    int resolution = (*data).size();
    int num_samples = (*data)[0].size();
    int dtype_size = (*data)[0][0].size();

    std::ofstream output_file(filename);

    // Write header
    output_file << resolution << "\t" << header << "\n";

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t";
        for (int n = 0; n < num_samples; n++) {
            for (int k = 0; k < dtype_size; k++) {
                output_file << (*data)[i][n][k] << "\t";
            }
        }
        output_file << "\n";
    }

    output_file.close();
}


template <class ModelType>
void generate_spin_configs(ModelType *model, std::vector<float> T, unsigned long long equilibration_steps, 
                                                              unsigned long long num_samples,
                                                              unsigned long long steps_per_sample,
                                                              unsigned long long num_threads, 
                                                              std::string foldername) {
    
    int resolution = T.size();

    ctpl::thread_pool threads(num_threads);
    std::vector<std::future<void>> results(resolution);

    auto spin_samples = [equilibration_steps, num_samples, steps_per_sample, model, foldername](int id, float T) {
        MonteCarlo<ModelType> *m = new MonteCarlo<ModelType>(model->clone());
        m->model->randomize_spins();

        m->steps(equilibration_steps, T);
        for (int j = 0; j < num_samples; j++) {
            m->steps(steps_per_sample, T);
            m->model->save_spins(foldername + "/Spins" + std::to_string(T) + "_" + std::to_string(j) + ".txt");
        }
    };

    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(spin_samples, T[i]);
    }

    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }
}   



