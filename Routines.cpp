#include "MonteCarlo.cpp"

#define BOLTZMANN_CONSTANT 0.08617

template <class ModelType>
void exchange(vector<MonteCarlo<ModelType>*> *models, vector<float> *T) {
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

template <class ModelType, class A>
vector<vector<A>> sample_pt(A sampling_func(ModelType*), ModelType *model, vector<float> *T, 
                                                                 unsigned long long steps_per_run, 
                                                                 unsigned long long num_samples, 
                                                                 unsigned long long steps_per_sample,
                                                                 int num_threads) {

    int resolution = T->size();

    vector<MonteCarlo<ModelType>*> models(resolution);

    for (int i = 0; i < resolution; i++) {
        models[i] = new MonteCarlo<ModelType>(model->clone());
        models[i]->model->randomize_spins();
    }

    vector<vector<A>> arr = vector<vector<A>>(resolution, vector<A>(num_samples));

    ctpl::thread_pool threads(num_threads);

    vector<future<void>> results(resolution);

    auto do_steps = [steps_per_run](int id, MonteCarlo<ModelType> *m, float T) {
        run_MC(m, steps_per_run, "trig", 3*T, T);
        //m->steps(steps_per_run, T);
    };

    // Do initial steps
    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(do_steps, models[i], (*T)[i]);
    }

    // Join threads
    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }

    auto take_samples = [steps_per_sample, sampling_func](int id, int i, int n,
                                                                  MonteCarlo<ModelType> *m, float T, 
                                                                  vector<A> *arr_i) {
        m->steps(steps_per_sample, T);
        A sample = sampling_func(m->model);
        (*arr_i)[n] = sample;
    };

    for (int n = 0; n < num_samples; n++) {
        for (int i = 0; i < resolution; i++) {
            results[i] = threads.push(take_samples, i, n, models[i], (*T)[i], &arr[i]);
        }

        // Join threads
        for (int i = 0; i < resolution; i++) {
            results[i].get();
        }

        exchange(&models, T);
    }

    return arr;
}

template <class ModelType, class A>
vector<vector<A>> sample_r(A sampling_func(ModelType*), ModelType *model, vector<float> *T, 
                                                        unsigned long long num_runs,
                                                        unsigned long long steps_per_run, 
                                                        unsigned long long num_samples, 
                                                        unsigned long long steps_per_sample,
                                                        int num_threads) {

    int resolution = T->size();

    vector<MonteCarlo<ModelType>*> models(resolution);

    for (int i = 0; i < resolution; i++) {
        models[i] = new MonteCarlo<ModelType>(model->clone());
        models[i]->model->randomize_spins();
    }

    vector<vector<A>> arr = vector<vector<A>>(resolution, vector<A>(num_runs*num_samples));

    ctpl::thread_pool threads(num_threads);

    vector<future<void>> results(resolution);

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
                                                                                    vector<A> *arr_i) {
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
vector<vector<vector<dtype>>> summary_statistics(vector<vector<vector<dtype>>> *data) {
    int resolution = (*data).size();
    int num_samples = (*data)[0].size();
    int dtype_size = (*data)[0][0].size();

    vector<vector<vector<dtype>>> arr(resolution, vector<vector<dtype>>(2, vector<dtype>(dtype_size)));

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
void write_data(vector<vector<vector<dtype>>> *data, vector<float> *T, string filename, string header = "") {
    int resolution = (*data).size();
    int dtype_size = (*data)[0][0].size();

    ofstream output_file(filename);

    // Write header
    output_file << resolution << "\t" << header << "\n";

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t";
        for (int k = 0; k < dtype_size; k++) {
            output_file << (*data)[i][0][k] << ", " << (*data)[i][1][k];
            if (k != dtype_size - 1) { output_file << "\t"; }
        }
        output_file << "\n";
    }

    output_file.close();
}

template <class ModelType>
void generate_spin_configs(ModelType *model, vector<float> T, unsigned long long equilibration_steps, 
                                                              unsigned long long num_samples,
                                                              unsigned long long steps_per_sample,
                                                              unsigned long long num_threads, 
                                                              string foldername) {
    
    int resolution = T.size();

    ctpl::thread_pool threads(num_threads);
    vector<future<void>> results(resolution);

    auto spin_samples = [equilibration_steps, num_samples, steps_per_sample, model, foldername](int id, float T) {
        MonteCarlo<ModelType> *m = new MonteCarlo<ModelType>(model->clone());
        m->model->randomize_spins();

        m->steps(equilibration_steps, T);
        for (int j = 0; j < num_samples; j++) {
            m->steps(steps_per_sample, T);
            m->model->save_spins(foldername + "/Spins" + to_string(T) + "_" + to_string(j) + ".txt");
        }
    };

    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(spin_samples, T[i]);
    }

    for (int i = 0; i < resolution; i++) {
        results[i].get();
    }
}   



