#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include<iostream>
#include<math.h>

using namespace std;
using namespace Eigen;


int main(int argc, char* argv[]) {
    
    string foldername = argv[1];

    // Load params config file
    ifstream paramfile(argv[2]);
    string Bs = argv[3];
    string line;

    getline(paramfile, line);
    getline(paramfile, line);
    vector<string> paramss = split(&line, ",");

    int N = stoi(paramss[0]);
    int L = stoi(paramss[1]);
    float J1 = stof(paramss[2]);
    float J2 = stof(paramss[3]);
    float K1 = stof(paramss[4]);
    float K2 = stof(paramss[5]);
    float K3 = stof(paramss[6]);

    getline(paramfile, line);
    getline(paramfile, line);
    getline(paramfile, line);
    paramss = split(&line, ",");

    int resolution = stoi(paramss[0]);
    int num_runs = stoi(paramss[1]);
    int steps_per_run = stoi(paramss[2]);
    int num_samples = stoi(paramss[3]);
    int steps_per_sample = stoi(paramss[4]);
    const int MCstep = N*N*L;



    Vector3f Bhat1; Bhat1 << 0.866025, 0.5, 0;
    Vector3f Bhat2; Bhat2 << 1., 0., 0.;
    Vector3f Bhat3; Bhat3 << 0., 0., 1.;
    float Bm = stof(Bs);
    Vector3f B1 = Bm*Bhat1;
    Vector3f B2 = Bm*Bhat2; 
    Vector3f B3 = Bm*Bhat3;
    
    float T_max = 30;
    float T_min = 0.2;
    float T;


    srand((unsigned)time( NULL ));
    TrigonalModel *model1 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, Bhat1);
    MonteCarlo<TrigonalModel> *m1 = new MonteCarlo(model1);
    TrigonalModel *model2 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, Bhat2);
    MonteCarlo<TrigonalModel> *m2 = new MonteCarlo(model2);
    TrigonalModel *model3 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, Bhat3);
    MonteCarlo<TrigonalModel> *m3 = new MonteCarlo(model3);

    Vector3f M1v; Vector3f M2v; Vector3f M3v;
    Vector3f M1v_avg; Vector3f M2v_avg; Vector3f M3v_avg;

    vector<float> M1(num_runs*num_samples); vector<float> M1squared(num_runs*num_samples);
    vector<float> M2(num_runs*num_samples); vector<float> M2squared(num_runs*num_samples);
    vector<float> M3(num_runs*num_samples); vector<float> M3squared(num_runs*num_samples);
    float M1_avg; float M2_avg; float M3_avg;
    float dM1; float dM2; float dM3;



    ofstream output_file1(foldername + "/SusceptibilityCurve1.txt");
    ofstream output_file2(foldername + "/SusceptibilityCurve2.txt");
    ofstream output_file3(foldername + "/SusceptibilityCurve3.txt");

    double seconds = 3.3*resolution*(num_runs*steps_per_run + num_samples*steps_per_sample)*N*N*L/3000000.;
    cout << "Expected completion time: " << seconds/60. << " minutes." << endl;

    auto start = chrono::high_resolution_clock::now();

    int ind;

    for (int i = 0; i < resolution; i++) {
        T = T_max*i/resolution + T_min*(resolution - i)/resolution;

        for (int j = 0; j < num_runs; j++) {
            model1->randomize_spins();
            model2->randomize_spins();
            model3->randomize_spins();
           
            run_MC(m1, steps_per_run*MCstep, "trig", J1, T);
            run_MC(m2, steps_per_run*MCstep, "trig", J1, T);
            run_MC(m3, steps_per_run*MCstep, "trig", J1, T);

            for (int k = 0; k < num_samples; k++) {
                run_MC(m1, steps_per_sample*MCstep, T);
                run_MC(m2, steps_per_sample*MCstep, T);
                run_MC(m3, steps_per_sample*MCstep, T);

                ind = j*num_samples + k;

                M1v = model1->get_magnetization();
                M1v_avg += M1v;
                M1[ind] = M1v.dot(Bhat1);
                M1squared[ind] = M1[ind]*M1[ind];

                M2v = model2->get_magnetization();
                M2v_avg += M2v;
                M2[ind] = M2v.dot(Bhat2);
                M2squared[ind] = M2[ind]*M2[ind];

                M3v = model3->get_magnetization();
                M3v_avg += M3v;
                M3[ind] = M3v.dot(Bhat3);
                M3squared[ind] = M3[ind]*M3[ind];
            }
        }

        M1_avg = avg(&M1);
        M2_avg = avg(&M2);
        M3_avg = avg(&M3);


        M1v_avg = M1v_avg/(num_runs*num_samples);
        M2v_avg = M2v_avg/(num_runs*num_samples);
        M3v_avg = M3v_avg/(num_runs*num_samples);

        dM1 = sqrt(abs(avg(&M1squared) - M1_avg*M1_avg));
        dM2 = sqrt(abs(avg(&M2squared) - M2_avg*M2_avg));
        dM3 = sqrt(abs(avg(&M3squared) - M3_avg*M3_avg));

        output_file1 << T << "\t" << "[" << M1v_avg[0] << "," << M1v_avg[1] << "," << M1v_avg[2] << "]\t" 
                     << M1_avg/Bm << "\t" << dM1/Bm << "\t" << model1->energy() << endl;
        output_file2 << T << "\t" << "[" << M2v_avg[0] << "," << M2v_avg[1] << "," << M2v_avg[2] << "]\t"
                     << M2_avg/Bm << "\t" << dM2/Bm << "\t" << model2->energy() << endl;
        output_file3 << T << "\t" << "[" << M3v_avg[0] << "," << M3v_avg[1] << "," << M3v_avg[2] << "]\t"
                     << M3_avg/Bm << "\t" << dM3/Bm << "\t" << model3->energy() << endl;
    }

    output_file1.close();
    output_file2.close();
    output_file3.close();

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    seconds = duration.count()/1000000.;

    cout << "Actual completion time: " << seconds/60. << " minutes." << endl;


}

