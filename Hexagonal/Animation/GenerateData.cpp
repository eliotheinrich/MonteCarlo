#include "../HexagonalModel.cpp"

int main(int argc, char* argv[]) {
    int N = stoi(argv[38]);
    float J1 = stof(argv[1]); float J2 = stof(argv[2]); float J3 = stof(argv[3]);
    float J4 = stof(argv[4]); float J5 = stof(argv[5]); float J6 = stof(argv[6]); 
    float J7 = stof(argv[7]); float J8 = stof(argv[8]); float J9 = stof(argv[9]);
    float D1 = stof(argv[10]); float D2 = stof(argv[11]); float D3 = stof(argv[12]); float phi = stof(argv[13]);
    float D4 = stof(argv[14]); float t4 = stof(argv[15]); float p4 = stof(argv[16]); 
    float D5 = stof(argv[17]); float t5 = stof(argv[18]); float p5 = stof(argv[19]);
    float D6 = stof(argv[20]); float t6 = stof(argv[21]); float p6 = stof(argv[22]); 
    float D7 = stof(argv[23]); float t7 = stof(argv[24]); float p7 = stof(argv[25]);
    float D8 = stof(argv[26]); float t8 = stof(argv[27]); float p8 = stof(argv[28]); 
    float D9 = stof(argv[29]); float t9 = stof(argv[30]); float p9 = stof(argv[31]);
    float K1 = stof(argv[32]); float K2 = stof(argv[33]); float K3 = stof(argv[34]);
    float B = stof(argv[35]);  float A = stof(argv[36]);  float S = stof(argv[37]);

    HexagonalModel model = HexagonalModel(N, J1, J2, J3, 
                                          J4, J5, J6, J7, J8, J9,
                                          D1, D2, D3, phi, 
                                          D4, t4, p4, D5, t5, p5, D6, t6, p6,
                                          D7, t7, p7, D8, t8, p8, D9, t9, p9,
                                          K1, K2, K3, 
                                          B, A, S);
    srand((unsigned)time( NULL ));

    int nsteps = stoi(argv[39]);
    float T_max = .05;
    float T_min = 0.02;

    run_MC<SpinTextureLogItem<2>>(&model, nsteps, "trig", T_max, T_min, true);

}
