#ifndef UTILITY_
#define UTILITY_

#include <iostream>
#include <cmath>
#include <chrono>
#include <random>

using namespace std;

int time_code(void fun()) {
    auto start = chrono::high_resolution_clock::now();

    fun();

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    return duration.count();
}

const inline int mod(int a, int b) {
  int c = a % b;
  return (c < 0) ? c + b : c;
}

template<typename T>
T avg(vector<T> *v) {
    int L = v->size();
    T F = T(0.);
    for (int i = 0; i < L; i++) {
        F = F + (*v)[i];
    }
    return F/float(L);
}

template<typename T>
T stdev(vector<T> *v, T av) {
    int L = v->size();
    T F = T(0.);
    T dF;
    for (int i = 0; i < L; i++) {
        dF = (*v)[i] - av;
        F = F + dF*dF;
    }
    return sqrt(F/float(L));
}

template<typename T>
T stdev(vector<T> *v) {
    T av = avg(v);
    return stdev<T>(v, av); 
}


vector<string> split(string *s, string delim) {
    vector<string> vals(0);

    string str = *s;
    int pos = 0;
    string token;
    while ((pos = str.find(delim)) != string::npos) {
        token = str.substr(0, pos);
        vals.push_back(token);
        str.erase(0, pos + delim.length());
    }

    vals.push_back(str);

    return vals;
}

class GaussianDist {
    private:
        minstd_rand rd;
        default_random_engine gen;
        normal_distribution<> dist;

    public:

        GaussianDist(float mean, float std) {
            this->rd.seed(rand());
            this->gen = default_random_engine(rd());
            this->dist = normal_distribution<>(mean, std);
        }

        float sample() {
            return dist(gen);
        }

};


 
#endif
