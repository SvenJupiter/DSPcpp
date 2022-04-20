#ifndef SJ_DSP_SIGNAL_HPP
#define SJ_DSP_SIGNAL_HPP

#include <vector> // 
#include <utility> // pair
#include "DSP/dsp_types.h" // real_t



namespace DSP {

namespace Discrete {

typedef std::vector<real_t> Signal;

// Returns the scalar dot product of u and v
DSP_FUNCTION real_t dot_product(const Signal& u, const Signal& v);

// Convolution and polynomial multiplication
// w = conv(u,v) returns the convolution of vectors u and v. 
// If u and v are vectors of polynomial coefficients, 
// convolving them is equivalent to multiplying the two polynomials.
DSP_FUNCTION Signal conv(const Signal& u, const Signal& v);


// Deconvolution and polynomial division
// [q,r] = deconv(u,v) deconvolves a vector v out of a vector u using long division, 
// and returns the quotient q and remainder r such that u = conv(v,q) + r. 
// If u and v are vectors of polynomial coefficients, then deconvolving them is equivalent to 
// dividing the polynomial represented by u by the polynomial represented by v.
DSP_FUNCTION std::pair<Signal, Signal> deconv(const Signal& u, const Signal& v);


// DSP_CLASS Signal: public std::vector<real_t> {

// }; // class Signal

}; // Namespace Discrete

}; // Namespace DSP

#endif // SJ_DSP_SIGNAL_HPP