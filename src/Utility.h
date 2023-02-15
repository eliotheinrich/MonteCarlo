#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <cmath>
#include <chrono>
#include <random>

int time_code(void fun());

const inline int mod(int a, int b) {
    int c = a % b;
    return (c < 0) ? c + b : c;
}

template<typename T>
T avg(std::vector<T> *v);

template<typename T>
T stdev(std::vector<T> *v, T av);

template<typename T>
T stdev(std::vector<T> *v);


std::vector<std::string> split(std::string *s, std::string delim);

class GaussianDist {
    private:
        std::minstd_rand rd;
        std::default_random_engine gen;
        std::normal_distribution<> dist;

    public:

        GaussianDist(float mean, float std);

        float sample();

};

 
#endif
