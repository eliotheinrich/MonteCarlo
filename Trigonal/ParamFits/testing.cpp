#include <iostream>
using namespace std;

class A {
    public:
        A() {
            cout << "A initialized!" << endl;
        }

        void foo() {
            cout << "Called foo from A" << endl;
        }
};

class B: public A {
    public:
        B() {
            cout << "B initialized!" << endl;
        }

        void foo() {
            cout << "Called foo from B" << endl;
            A::foo();
        }
};

int main() {
    //A *a = new A();
    //a->foo();

    B *b = new B();
    b->foo();
}
