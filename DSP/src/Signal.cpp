#include "DSP/Discrete/Signal.hpp"

using namespace DSP::Discrete;

real_t dot_product(const Signal& u, const Signal& v) {
    if (u.size() != v.size() ) { return 0;}
    const size_t size = u.size();

    // Init
    real_t sum = v[0] * u[0];

    // Add
    for (size_t k = 1; k < size; ++k) {
        sum += v[k] * u[k];
    }

    // Return
    return sum;
}


Signal conv(const Signal& u, const Signal& v) {
    //Check
    if (u.size() == 0 || v.size() == 0) { return Signal(); }
    const size_t conv_size = (u.size() + v.size() - 1);

    // From Matlab:
    // w = conv(u,v) returns the convolution of vectors u and v. 
    // If u and v are vectors of polynomial coefficients, 
    // convolving them is equivalent to multiplying the two polynomials.
    size_t start, end;
    Signal w(conv_size, 0);
    for (size_t	k = 0; k < conv_size; k++) {

        // Determine start and end
        start = (k < u.size() ? 0 : k - (u.size() - 1));
        end = (k < v.size() ? k : (v.size() - 1));

        // // Init
        // w[k] = 0;

        // Calculate k-th element
        for (size_t j = start; j <= end; ++j) {
            w[k] += v[j] * u[k-j];
        }
    }

    return w;
}


std::pair<Signal, Signal> deconv(const Signal& u, const Signal& v) {

    // Check
    if (v[0] == 0) { return std::make_pair(Signal(), Signal()); }
    if (u.size() < v.size()) { return std::make_pair(Signal(), Signal());; }        

    // calculate sizes
    const size_t q_size = u.size() - v.size() + 1;
    // const size_t r_size = v.size() - 1;

    // we copy 'u' into 'r' at the start
    Signal r = u;
    Signal q(q_size, 0);

    // Deconvolution 
    for (size_t k = 0; k < q_size; ++k) {

        // Divide
        q[k] = r[k] / v[0];

        // Subtract
        r[k] = 0;
        for (size_t j = 1; j < v.size(); ++j) {
            r[k+j] -= q[k] * v[j]; 
        }
    }

    return std::make_pair(q, r);
}