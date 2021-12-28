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

class LatticeIterator {
    public:
        int n1; int n2; int n3; int s;

        int N1; int N2; int N3; int sl;

        int counter;
        int t;

        LatticeIterator(int N1, int N2, int N3, int sl) {
            this->n1 = 0; this->n2 = 0; this->n3 = 0; this->s = 0;

            this->N1 = N1; this->N2 = N2; this->N3 = N3; this->sl = sl;

            this->counter = 0;
        }

        void next() {
            counter++;

            t = counter;
            s = t % sl;
            t = (t - s)/sl;
            n1 = t % N1;
            t = (t - n1)/N1;
            n2 = t % N2;
            t = (t - n2)/N2;
            n3 = t % N3;

            if (counter == N1*N2*N3*sl) {
                counter = 1;
            }
        }

        void reset() {
            n1 = 0; n2 = 0; n3 = 0; sl = 0;
            counter = 0;
        }
};

 
#endif
