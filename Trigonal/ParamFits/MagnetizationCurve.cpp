#include "../TrigonalModel.cpp"
#include "../../Utility.cpp"
#include<iostream>
#include<math.h>

using namespace std;
using namespace Eigen;


int main(int argc, char* argv[]) {
    
    string foldername = argv[1];
    string Ts = argv[3];

    // Load params config file
    ifstream paramfile(argv[2]);
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



    float T = stof(Ts);

    Vector3f Bhat1; Bhat1 << 0., 1., 0;
    Vector3f Bhat2; Bhat2 << 1., 0., 0.;
    Vector3f B1; Vector3f B2;
    float Bm;
    


    srand((unsigned)time( NULL ));
    TrigonalModel *model1 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, Bhat1);
    MonteCarlo<TrigonalModel> *m1 = new MonteCarlo(model1);
    TrigonalModel *model2 = new TrigonalModel(N, L, J1, J2, K1, K2, K3, Bhat2);
    MonteCarlo<TrigonalModel> *m2 = new MonteCarlo(model2);

    Vector3f M1v; Vector3f M2v;
    Vector3f M1v_avg; Vector3f M2v_avg;

    vector<float> M1(num_runs*num_samples); vector<float> M1squared(num_runs*num_samples);
    vector<float> M2(num_runs*num_samples); vector<float> M2squared(num_runs*num_samples);
    float M1_avg; float M2_avg;
    float dM1; float dM2;

    const float B_max = 1.5;
    const float B_min = 0.;



    ofstream output_file1(foldername + "/MagnetizationCurve1.txt");
    ofstream output_file2(foldername + "/MagnetizationCurve2.txt");

    double seconds = 2.*resolution*num_runs*(steps_per_run + num_samples*steps_per_sample)*N*N*L/3000000.;
    cout << "Expected completion time: " << seconds/60. << " minutes." << endl;

    auto start = chrono::high_resolution_clock::now();

    int ind;

    for (int i = 0; i < resolution; i++) {
        Bm = B_max*i/resolution + B_min*(resolution - i)/resolution;
        B1 = Bm*Bhat1;
        B2 = Bm*Bhat2;
        model1->B = B1;
        model2->B = B2;

        for (int j = 0; j < num_runs; j++) {
            model1->randomize_spins();
            model2->randomize_spins();
           
            run_MC(m1, steps_per_run*MCstep, "trig", J1, T);
            run_MC(m2, steps_per_run*MCstep, "trig", J1, T);

            for (int k = 0; k < num_samples; k++) {
                run_MC(m1, steps_per_sample*MCstep, T);
                run_MC(m2, steps_per_sample*MCstep, T);

                ind = j*num_samples + k;

                M1v = model1->get_magnetization();
                M1v_avg += M1v;
                M1[ind] = model1->get_magnetization().dot(Bhat1);
                M1squared[ind] = M1[ind]*M1[ind];

                M2v = model2->get_magnetization();
                M2v_avg += M2v;
                M2[ind] = model2->get_magnetization().dot(Bhat2);
                M2squared[ind] = M2[ind]*M2[ind];
            }
        }

//        model1->save_spins(foldername + "/ab" + to_string(Bm) + ".txt");
//        model2->save_spins(foldername + "/a" + to_string(Bm) + ".txt");

        M1_avg = avg(&M1);
        M2_avg = avg(&M2);

        M1v_avg = M1v_avg/(num_runs*num_samples);
        M2v_avg = M2v_avg/(num_runs*num_samples);

        dM1 = sqrt(abs(avg(&M1squared) - M1_avg*M1_avg));
        dM2 = sqrt(abs(avg(&M2squared) - M2_avg*M2_avg));

        output_file1 << Bm << "\t" << "[" << M1v_avg[0] << "," << M1v_avg[1] << "," << M1v_avg[2] << "]\t" 
                     << M1_avg << "\t" << dM1 << "\t" << model1->energy() << endl;
        output_file2 << Bm << "\t" << "[" << M2v_avg[0] << "," << M2v_avg[1] << "," << M2v_avg[2] << "]\t"
                     << M2_avg << "\t" << dM2 << "\t" << model2->energy() << endl;
    }

    output_file1.close();
    output_file2.close();

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    seconds = duration.count()/1000000.;

    cout << "Actual completion time: " << seconds/60. << " minutes." << endl;


}

