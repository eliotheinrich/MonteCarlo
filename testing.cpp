#include <thread>
#include <ctpl_stl.h>
#include <iostream>
#include <vector>

using namespace std;

int func(int id, float f) {
    cout << "id = " << id << ", f = " << f << endl;
    int s = 0;
    for (int i = 0; i < 1000; i++) {
        s += i;
    }
    return s;
}

int main() {
    ctpl::thread_pool pool(4);

    int num_jobs = 20;
    vector<future<int>> results(num_jobs);

    for (int i = 0; i < num_jobs; i++) {
        results[i] = pool.push(func, 0.1);
    }

    for (int i = 0; i < num_jobs; i++) {
        cout << results[i].get() << endl;
    }
}
