#include "../HexagonalModel.cpp"

int main() {
    int N = 20;
    float JNN = 1.;
    float JNNN = 1.;
    float DNN = 0.2;
    float DNNN = 0.1;
    float p = 0.;
    float t = 0.;
    float KNN = 1.2;

    float J1 = JNN; float J2 = JNN; float J3 = JNN;
    float J4 = JNNN*1.1; float J5 = JNNN*2; float J6 = JNNN*6; float J7 = JNNN*0.5; float J8 = JNNN*2.5; float J9 = JNNN*0.1;
    float D1 = DNN; float D2 = DNN; float D3 = DNN; float phi = 0.;
    float D4 = DNNN; float D5 = DNNN; float D6 = DNNN; float D7 = DNNN; float D8 = DNNN; float D9 = DNNN;
    float t4 = t; float t5 = t; float t6 = t; float t7 = t; float t8 = t; float t9 = t;
    float p4 = p; float p5 = p; float p6 = p; float p7 = p; float p8 = p; float p9 = p;
    float K1 = KNN; float K2 = KNN; float K3 = KNN;
    float B = 0.25;  float A = 0.3;  float S = 1.;
    HexagonalModel model = HexagonalModel(N, J1, J2, J3, 
                                          J4, J5, J6, J7, J8, J9,
                                          D1, D2, D3, phi, 
                                          D4, t4, p4, D5, t5, p5, D6, t6, p6,
                                          D7, t7, p7, D8, t8, p8, D9, t9, p9,
                                          K1, K2, K3, 
                                          B, A, S);

    MonteCarlo<HexagonalModel, SpinTextureLogItem<2>> m(&model);

    int nsteps = 10000000;
    int num_log_updates = 200;
    int num_temp_updates = 100;
    string cooling_schedule = "trig";
    float T_max = .2;
    float T_min = 0.0001;
    string filename = "Log.txt";
    m.steps(nsteps, T_max, T_min, cooling_schedule, filename, num_temp_updates, num_log_updates);
    cout << model.get_magnetization() << endl;
//    m.steps(nsteps, T_min, T_min, "const", filename, num_temp_updates, num_log_updates);
}
