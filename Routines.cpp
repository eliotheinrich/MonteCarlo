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

template <class ModelType, class dtype>
void sample_pt(vector<dtype> sampling_func(ModelType*), ModelType *model, vector<float> *T, 
                                                        int steps_per_run, 
                                                        int num_samples, 
                                                        int steps_per_sample,
                                                        int num_threads,
                                                        string filename) {

    int resolution = T->size();

    vector<MonteCarlo<ModelType>*> models(resolution);

    for (int i = 0; i < resolution; i++) {
        models[i] = new MonteCarlo<ModelType>(model->clone());
        models[i]->model->randomize_spins();
    }

    int dtype_size = sampling_func(model).size();
    vector<vector<vector<dtype>>> arr = vector<vector<vector<dtype>>>(resolution, vector<vector<dtype>>(dtype_size, vector<dtype>(num_samples)));

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

    auto take_samples = [steps_per_sample, dtype_size, sampling_func](int id, int i, int n,
                                                                         MonteCarlo<ModelType> *m, float T, 
                                                                         vector<vector<dtype>> *arr_i) {
        m->steps(steps_per_sample, T);
        auto sample = sampling_func(m->model);
        for (int k = 0; k < dtype_size; k++) {
            (*arr_i)[k][n] = sample[k];
        }
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

    vector<vector<dtype>> avgs(resolution, vector<dtype>(dtype_size));
    vector<vector<dtype>> errs(resolution, vector<dtype>(dtype_size));

    for (int i = 0; i < resolution; i++) {
        for (int k = 0; k < dtype_size; k++) {
            avgs[i][k] = avg(&arr[i][k]);
            errs[i][k] = stdev(&arr[i][k], avgs[i][k]);
        }
    }

    ofstream output_file(filename);

    // Write header
    output_file << resolution << endl;

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t";
        for (int k = 0; k < dtype_size; k++) {
            output_file << avgs[i][k] << ", " << errs[i][k];
            if (k != dtype_size - 1) { output_file << "\t"; }
        }
        output_file << "\n";
    }

    output_file.close();
}

template <class ModelType, class dtype>
void sample_r(vector<dtype> sampling_func(ModelType*), ModelType *model, vector<float> *T, 
                                                        int num_runs,
                                                        int steps_per_run, 
                                                        int num_samples, 
                                                        int steps_per_sample,
                                                        int num_threads,
                                                        string filename) {

    int resolution = T->size();

    vector<MonteCarlo<ModelType>*> models(resolution);

    for (int i = 0; i < resolution; i++) {
        models[i] = new MonteCarlo<ModelType>(model->clone());
        models[i]->model->randomize_spins();
    }

    int dtype_size = sampling_func(model).size();
    vector<vector<vector<dtype>>> arr = vector<vector<vector<dtype>>>(resolution, 
                                               vector<vector<dtype>>(dtype_size, 
                                                      vector<dtype>(num_runs*num_samples)));

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

    auto take_samples = [num_samples, steps_per_sample, num_runs, steps_per_run, dtype_size, sampling_func](int id, int i,
                                                                                   MonteCarlo<ModelType> *m, float T, 
                                                                                   vector<vector<dtype>> *arr_i) {
        for (int n = 0; n < num_runs; n++) {
            m->model->randomize_spins();
            m->steps(steps_per_run, T);
            for (int j = 0; j < num_samples; j++) {
                auto sample = sampling_func(m->model);
                for (int k = 0; k < dtype_size; k++) {
                    (*arr_i)[k][n*num_samples + j] = sample[k];
                }
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

    vector<vector<dtype>> avgs(resolution, vector<dtype>(dtype_size));
    vector<vector<dtype>> errs(resolution, vector<dtype>(dtype_size));

    for (int i = 0; i < resolution; i++) {
        for (int k = 0; k < dtype_size; k++) {
            avgs[i][k] = avg(&arr[i][k]);
            errs[i][k] = stdev(&arr[i][k], avgs[i][k]);
        }
    }

    ofstream output_file(filename);

    // Write header
    output_file << resolution << endl;

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t";
        for (int k = 0; k < dtype_size; k++) {
            output_file << avgs[i][k] << ", " << errs[i][k];
            if (k != dtype_size - 1) { output_file << "\t"; }
        }
        output_file << "\n";
    }

    output_file.close();
}

template <class ModelType>
void generate_spin_configs(ModelType *model, vector<float> T, int equilibration_steps, 
                                                              int num_samples,
                                                              int steps_per_sample,
                                                              int num_threads, 
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



