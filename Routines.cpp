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


template <class ModelType>
void stiffness_run_pt(ModelType *model, vector<float> *T, int steps_per_run, 
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

    vector<vector<float>> dE = vector<vector<float>>(resolution, vector<float>(num_samples));
    vector<vector<float>> ddE = vector<vector<float>>(resolution, vector<float>(num_samples));

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

    auto twist_samples = [steps_per_sample](int id, int i, int n, 
                                            MonteCarlo<ModelType> *m, float T, 
                                            vector<float> *dE_samples, vector<float> *ddE_samples) {
        m->steps(steps_per_sample, T);
        auto sample = m->model->twist_stiffness();
        (*dE_samples)[n] = sample[0];
        (*ddE_samples)[n] = sample[1];
    };

    for (int n = 0; n < num_samples; n++) {
        for (int i = 0; i < resolution; i++) {
            results[i] = threads.push(twist_samples, i, n, models[i], (*T)[i], &dE[i], &ddE[i]);
        }

        // Join threads
        for (int i = 0; i < resolution; i++) {
            results[i].get();
        }

        exchange(&models, T);
    }

    vector<float> avg_dE(resolution);
    vector<float> avg_ddE(resolution);

    vector<float> err_dE(resolution);
    vector<float> err_ddE(resolution);

    for (int i = 0; i < resolution; i++) {
        avg_dE[i] = avg(&dE[i]);
        avg_ddE[i] = avg(&ddE[i]);
        err_dE[i] = stdev(&dE[i]);
        err_ddE[i] = stdev(&ddE[i]);
    }

    ofstream output_file(filename);

    // Write header
    output_file << resolution << endl;

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t" << avg_dE[i] << "\t" << err_dE[i] << "\t" << avg_ddE[i] << "\t" << err_ddE[i] << endl;
    }

    output_file.close();
}

template <class ModelType>
void stiffness_run(ModelType *model, vector<float> *T, int num_runs,
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

    vector<vector<float>> dE = vector<vector<float>>(resolution, vector<float>(num_samples*num_runs));
    vector<vector<float>> ddE = vector<vector<float>>(resolution, vector<float>(num_samples*num_runs));

    ctpl::thread_pool threads(num_threads);

    vector<future<void>> results(resolution);

    auto twist_samples = [num_runs, steps_per_run, num_samples, steps_per_sample](int id, int i, 
                                                                                  MonteCarlo<ModelType> *m, float T, 
                                                                                  vector<float> *dE_samples, 
                                                                                  vector<float> *ddE_samples) {
        for (int n = 0; n < num_runs; n++) {
            m->model->randomize_spins();
            m->steps(steps_per_run, T);
            for (int i = 0; i < num_samples; i++) {
                auto sample = m->model->twist_stiffness();
                m->steps(steps_per_sample, T);
                (*dE_samples)[n] = sample[0];
                (*ddE_samples)[n] = sample[1];
            }
        }
    };

    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(twist_samples, i, models[i], (*T)[i], &dE[i], &ddE[i]);
    }

    vector<float> avg_dE(resolution);
    vector<float> avg_ddE(resolution);

    vector<float> err_dE(resolution);
    vector<float> err_ddE(resolution);

    for (int i = 0; i < resolution; i++) {
        avg_dE[i] = avg(&dE[i]);
        avg_ddE[i] = avg(&ddE[i]);
        err_dE[i] = stdev(&dE[i]);
        err_ddE[i] = stdev(&ddE[i]);
    }

    ofstream output_file(filename);

    // Write header
    output_file << resolution << endl;

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t" << avg_dE[i] << "\t" << err_dE[i] << "\t" << avg_ddE[i] << "\t" << err_ddE[i] << endl;
    }

    output_file.close();
}

template <class ModelType>
void susceptibility_run_pt(ModelType *model, vector<float> *T, int steps_per_run, 
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

    vector<vector<float>> X = vector<vector<float>>(resolution, vector<float>(num_samples));
    vector<vector<float>> E = vector<vector<float>>(resolution, vector<float>(num_samples));

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

    auto susceptibility_samples = [steps_per_sample](int id, int i, int n, 
                                                     MonteCarlo<ModelType> *m, float T, 
                                                     vector<float> *X, vector<float> *E) {
        m->steps(steps_per_sample, T);
        (*X)[n] = m->model->get_magnetization().norm()/pow(m->model->B.norm(), 2);
        (*E)[n] = m->model->energy()/(m->model->N*m->model->N*m->model->L);
    };

    for (int n = 0; n < num_samples; n++) {
        for (int i = 0; i < resolution; i++) {
            results[i] = threads.push(susceptibility_samples, i, n, models[i], (*T)[i], &X[i], &E[i]);
        }

        // Join threads
        for (int i = 0; i < resolution; i++) {
            results[i].get();
        }

        exchange(&models, T);
    }
    

    vector<float> avg_X(resolution);
    vector<float> avg_E(resolution);

    vector<float> err_X(resolution);
    vector<float> err_E(resolution);

    for (int i = 0; i < resolution; i++) {
        avg_X[i] = avg(&X[i]);
        avg_E[i] = avg(&E[i]);
        err_X[i] = stdev(&X[i], avg_X[i]);
        err_E[i] = stdev(&E[i], avg_E[i]);
    }

    ofstream output_file(filename);

    // Write header
    output_file << resolution << endl;

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t" << avg_X[i] << "\t" << err_X[i] << "\t" << avg_E[i] << "\t" << err_E[i] << endl;
    }

    output_file.close();
}

template <class ModelType>
void susceptibility_run(ModelType *model, vector<float> *T, int num_runs,
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

    vector<vector<float>> X = vector<vector<float>>(resolution, vector<float>(num_samples*num_runs));
    vector<vector<float>> E = vector<vector<float>>(resolution, vector<float>(num_samples*num_runs));

    ctpl::thread_pool threads(num_threads);

    vector<future<void>> results(resolution);

    auto susceptibility_samples = [num_runs, steps_per_run, num_samples, steps_per_sample](int id, int i, 
                                                                                           MonteCarlo<ModelType> *m, float T, 
                                                                                           vector<float> *X, 
                                                                                           vector<float> *E) {
        for (int n = 0; n < num_runs; n++) {
            m->model->randomize_spins();
            m->steps(steps_per_run, T);
            for (int i = 0; i < num_samples; i++) {
                (*X)[n*num_samples + i] = m->model->get_magnetization().norm()/pow(m->model->B.norm(), 2);
                (*E)[n*num_samples + i] = m->model->energy()/(m->model->V);
                m->steps(steps_per_sample, T);
            }
        }
    };

    for (int i = 0; i < resolution; i++) {
        results[i] = threads.push(susceptibility_samples, i, models[i], (*T)[i], &X[i], &E[i]);
    }

    vector<float> avg_X(resolution);
    vector<float> avg_E(resolution);

    vector<float> err_X(resolution);
    vector<float> err_E(resolution);

    for (int i = 0; i < resolution; i++) {
        avg_X[i] = avg(&X[i]);
        avg_E[i] = avg(&E[i]);
        err_X[i] = stdev(&X[i], avg_X[i]);
        err_E[i] = stdev(&E[i], avg_E[i]);
    }

    ofstream output_file(filename);

    // Write header
    output_file << resolution << endl;

    for (int i = 0; i < resolution; i++) {
        output_file << (*T)[i] << "\t" << avg_X[i] << "\t" << err_X[i] << "\t" << avg_E[i] << "\t" << err_E[i] << endl;
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



