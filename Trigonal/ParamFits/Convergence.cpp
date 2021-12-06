#include "../TrigonalModel.cpp"
#include<iostream>

using namespace std;
using namespace Eigen;


int main() {
    // Load params config file
    ifstream paramfile("params.txt");
    string line;

    getline(paramfile, line);
    getline(paramfile, line);
    vector<string> paramss = split(&line, ",");

    int N = stoi(paramss[0]);
    int L = stoi(paramss[1]);
    float J1 = stof(paramss[2]);
    float J2 = stof(paramss[3]);
    float A = stof(paramss[4]);
    float K = stof(paramss[5]);
    float S = stof(paramss[6]);


    const float T = 0.155;

    Vector3f Bhat1; Bhat1 << 0.866025, 0.5, 0;
    Vector3f Bhat2; Bhat2 << 1., 0., 0.;
    Vector3f B1; Vector3f B2;
    float Bm;
    

    const int MCstep = N*N*L;

    srand((unsigned)time( NULL ));
    TrigonalModel *model1;
    TrigonalModel *model2;

    model1 = new TrigonalModel(N, L, J1, J2, A, K, 0.1*Bhat1, S);
    model2 = new TrigonalModel(N, L, J1, J2, A, K, 0.1*Bhat2, S);

//    cout << "Model 1 took " << equilibrate(model1, T)/MCstep << " MCsteps to equilibrate." << endl;
//    cout << "Model 2 took " << equilibrate(model2, T)/MCstep << " MCsteps to equilibrate." << endl;
//

    run_MC<EnergyLogItem>(model1, 3000*MCstep, "trig", J1, T, true, 100, "log1.txt");
    run_MC<EnergyLogItem>(model2, 3000*MCstep, "trig", J1, T, true, 100, "log2.txt");

    cout << "Min1: " << model1->energy() << endl;
    cout << "Min2: " << model2->energy() << endl;

}

