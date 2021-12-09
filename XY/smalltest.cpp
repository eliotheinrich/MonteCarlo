#include <iostream>
#include <chrono>
#include <vector>
#include <thread>

using namespace std;

int get_int(int i) {
    return int(float(i));
}


int main() {
    int num_threads = 10000;
    vector<thread> threads(num_threads);
    vector<int> ints(num_threads);

    for (int i = 0; i < num_threads; i++) {
        threads[i] = thread(get_int, i);
    }

    for (int i = 0; i < num_threads; i++) {
        threads[i].join();
    }

    cout << "Finished!" << endl;
}
