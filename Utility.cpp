#ifndef UTILITY_
#define UTILITY_

#include <iostream>
#include <cmath>
#include <chrono>
#include <random>

int time_code(void fun()) {
    auto start = std::chrono::high_resolution_clock::now();

    fun();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    return duration.count();
}

const inline int mod(int a, int b) {
  int c = a % b;
  return (c < 0) ? c + b : c;
}

template<typename T>
T avg(std::vector<T> *v) {
    int L = v->size();
    T F = T(0.);
    for (int i = 0; i < L; i++) {
        F = F + (*v)[i];
    }
    return F/float(L);
}

template<typename T>
T stdev(std::vector<T> *v, T av) {
    int L = v->size();
    T F = T(0.);
    T dF;
    for (int i = 0; i < L; i++) {
        dF = (*v)[i] - av;
        F = F + dF*dF;
    }
    return std::sqrt(F/float(L));
}

template<typename T>
T stdev(std::vector<T> *v) {
    T av = avg(v);
    return stdev<T>(v, av); 
}


std::vector<std::string> split(std::string *s, std::string delim) {
    std::vector<std::string> vals(0);

    std::string str = *s;
    int pos = 0;
    std::string token;
    while ((pos = str.find(delim)) != std::string::npos) {
        token = str.substr(0, pos);
        vals.push_back(token);
        str.erase(0, pos + delim.length());
    }

    vals.push_back(str);

    return vals;
}

class GaussianDist {
    private:
        std::minstd_rand rd;
        std::default_random_engine gen;
        std::normal_distribution<> dist;

    public:

        GaussianDist(float mean, float std) {
            this->rd.seed(rand());
            this->gen = std::default_random_engine(rd());
            this->dist = std::normal_distribution<>(mean, std);
        }

        float sample() {
            return dist(gen);
        }

};

 
#endif
