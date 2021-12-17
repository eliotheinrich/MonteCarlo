#include <iostream>
#include <random>

using namespace std;


int main() {
    minstd_rand r;
    r.seed(10);
    for (int i = 0; i < 10000; i++) {
        cout << float(r())/RAND_MAX << endl;
    }

}
