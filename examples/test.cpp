#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "DSP/Math/Polynomial.hpp"
#include "DSP/Math/Matrix.hpp"
#include "DSP/Discrete/Signal.hpp"
#include "DSP/Discrete/pidController.hpp"
using namespace std;
using namespace DSP::Math;
using namespace DSP::Discrete;





void export_to_csv(const string& filename, const Signal& t, const Signal& u, const Signal& y) {

    // Open file
    ofstream file;
    file.open(filename);

    // Write header
    file << "t, u, y" << endl;

    // Write data
    for (size_t k = 0; k < t.size(); ++k) {
        file << t[k] << ", " << u[k] << ", " << y[k] << endl;
    }

    // Close file
    file.close();
}

void export_to_csv(const string& filename, const Signal& t, const Signal& w, const Signal& e, const Signal& y, const Signal& x) {

    // Open file
    ofstream file;
    file.open(filename);

    // Write header
    file << "t, w, e, y, x" << endl;

    // Write data
    for (size_t k = 0; k < t.size(); ++k) {
        file << t[k] << ", " << w[k] << ", " << e[k] << ", " << y[k] << ", " << x[k] << endl;
    }

    // Close file
    file.close();
}

real_t step(const real_t t, const real_t StepTime = 0) {
    return ((t >= StepTime) ? 1 : 0);
}

void test_polymomial_1() {

    Polynomial p = {-4, 6, -3, 1};
    cout << "p = " << p.to_string() << endl;

    cout << "P: ";
    for (auto& ak : p) {
        ak = 42;
    }
    cout << endl;

    cout << "p = " << p.to_string() << endl;

}

void test_polymomial_2() {  

    const vector<real_t> x = {1, 2, -1, 0, 5};
    const vector<real_t> y = {0.2, 0.4, -3, -0.5, 0.6};

    const auto V = Matrix::vandermonde(2, x);
    cout << "V: " << endl << V.to_string() << endl;
    

    const auto p = Polynomial::fit(2, x, y);
    cout << "p: " << p.to_string() << endl;

}

void test_matrix() {

    const real_t Amat[] = { 
        1, 4, 2,
        -3, -5, -2,
        1, -1, 0,
        1, 7, 4,
        -9, -4, 1
    };

    const real_t bvec[] = {1, -2, 4, 5, -2};

    const Matrix A(5, 3, Amat);
    const Matrix b(5, 1, bvec);
    const Matrix x = Matrix::solve(A, b);

    cout << "A: " << endl << A.to_string() << endl;
    cout << "~A: " << endl << (~A).to_string() << endl;

    cout << "b: " << endl << b.to_string() << endl;
    cout << "x: " << endl << x.to_string() << endl;
    cout << "(b / A): " << endl << (b / A).to_string() << endl;
    cout << "pinv(A) * b: " << endl << (A.pinv() * b).to_string() << endl;
    
}

void test_zss_pt1() {

    const real_t Tstep = 0.0;
    const real_t Ts = 0.01;
    const real_t duration = 30;
    // const size_t num_samples = 1 + round(duration / Ts);
    
    // PT1
    auto pt1 = zStateSpace::PT1(2, 3, Ts);


    // Simulate
    real_t tk, uk, yk;
    Signal t, u, y;
    for (tk = 0; tk <= duration; tk += Ts) {
        uk = step(tk);
        yk = pt1(uk);
        t.push_back(tk);
        u.push_back(uk);
        y.push_back(yk);
    }

    // Export
    export_to_csv("PT1-Test.csv", t, u, y);

}

void test_zss_hpf() {

    const real_t Tstep = 0.0;
    const real_t Ts = 0.1;
    const real_t duration = 10;
    // const size_t num_samples = 1 + round(duration / Ts);
    
    // Highpass-Filter
    auto hpf = zStateSpace::HighpassFilter(1, 1, Ts);


    // Simulate
    real_t tk, uk, yk;
    Signal t, u, y;
    for (tk = 0; tk <= duration; tk += Ts) {
        uk = step(tk);
        yk = hpf(uk);
        t.push_back(tk);
        u.push_back(uk);
        y.push_back(yk);
    }

    // Export
    export_to_csv("HPF-Test.csv", t, u, y);

}


void test_pid() {

    const real_t Tstep = 0.0;
    const real_t Ts = 0.1;
    const real_t duration = 30;
    // const size_t num_samples = 1 + round(duration / Ts);
    
    // PT1
    auto pt1 = zStateSpace::PT1(2, 3, Ts);

    // PID
    auto pid = pidController(1, 1, 0, 0.1, Ts);

    // Simulate
    real_t tk, wk, ek, yk, xk = 0;
    Signal t, w, e, y, x;
    for (tk = 0; tk <= duration; tk += Ts) {
        wk = step(tk);
        ek = wk - xk; // e[k] = w[k] - x[k-1]
        yk = pid(ek);
        xk = pt1(yk);
        t.push_back(tk);
        w.push_back(wk);
        e.push_back(ek);
        y.push_back(yk);
        x.push_back(xk);
    }

    // Export
    export_to_csv("PID-Test.csv", t, w, e, y, x);

}


void test_pid2() {

    const real_t Tstep = 0.0;
    const real_t Ts = 0.1;
    const real_t duration = 20;
    // const size_t num_samples = 1 + round(duration / Ts);
    
    // Highpass-Filter
    auto hpf = zStateSpace::HighpassFilter(1, 2, Ts);


    // PID
    auto pid = pidController(0, 1, 0, 0.1, Ts);

    // Simulate
    real_t tk, wk, ek, yk, xk = 0;
    Signal t, w, e, y, x;
    for (tk = 0; tk <= duration; tk += Ts) {
        wk = step(tk);
        ek = hpf(wk);
        yk = pid(ek);
        xk = yk;
        t.push_back(tk);
        w.push_back(wk);
        e.push_back(ek);
        y.push_back(yk);
        x.push_back(xk);
    }

    // Export
    export_to_csv("PID-Test_v2.csv", t, w, e, y, x);

}


int main() {

    cout << "Hello World!" << endl;
    // test_polymomial_1();
    // test_polymomial_2();
    // test_matrix();
    // test_zss_pt1();
    // test_zss_hpf();
    // test_pid();
    test_pid2();
    
    cout << "Bye bye..." << endl;
}
