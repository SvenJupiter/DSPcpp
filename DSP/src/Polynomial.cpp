#include <cstdlib> // malloc, realloc, free
#include <string> // memcpy, memset, memmove
#include <sstream> // stringstream

#include "DSP/Math/Polynomial.hpp"
#include "DSP/Math/Matrix.hpp"

using namespace DSP::Math;

// Constructor
Polynomial::Polynomial():
order(0),
coeffs(nullptr)
{

}

// Constructor
Polynomial::Polynomial(const size_t order):
order(0),
coeffs(nullptr)
{
    this->coeffs = (real_t*) malloc((order + 1) * sizeof(real_t));
    if (this->coeffs != nullptr) {
        this->order = order;
    }
}

// Constructor
Polynomial::Polynomial(const size_t other_order, const real_t* const other_coeffs):
order(0),
coeffs(nullptr)
{
    this->coeffs = (real_t*) malloc((other_order + 1) * sizeof(real_t));
    if (this->coeffs != nullptr) {
        this->order = other_order;
        if (other_coeffs != nullptr) {
            memcpy(this->coeffs, other_coeffs, (other_order + 1) * sizeof(real_t));
        }
        else {
            memset(this->coeffs, 0, (other_order + 1) * sizeof(real_t));
        }
    }
}

// Constructor
Polynomial::Polynomial(const std::vector<real_t>& other_coeffs):
order(0),
coeffs(nullptr)
{
    if (other_coeffs.size() > 0) {
        this->coeffs = (real_t*) malloc(other_coeffs.size() * sizeof(real_t));
        if (this->coeffs != nullptr) {
            this->order = other_coeffs.size() - 1;
            memcpy(this->coeffs, other_coeffs.data(), other_coeffs.size() * sizeof(real_t));
        }
    }
}

// Initilize
Polynomial::Polynomial(const std::initializer_list<real_t>& il):
order(0),
coeffs(nullptr)
{
    this->coeffs = (real_t*) malloc(il.size() * sizeof(real_t));
    if (this->coeffs != nullptr) {
        this->order = il.size() - 1;
        size_t k = 0;
        for (const auto& element : il) {
            this->coeffs[k] = element;
            ++k;
        }
    }
}

// Initilize
void Polynomial::operator=(const std::initializer_list<real_t>& il) {
    if (this->coeffs != nullptr) { free(this->coeffs); this->order = 0; this->coeffs = nullptr; }
    this->coeffs = (real_t*) malloc(il.size() * sizeof(real_t));
    if (this->coeffs != nullptr) {
        this->order = il.size() - 1;
        size_t k = 0;
        for (const auto& element : il) {
            this->coeffs[k] = element;
            ++k;
        }
    }
}


// copy
Polynomial::Polynomial(const Polynomial& other):
order(0),
coeffs(nullptr)
{
    if (other.coeffs != nullptr) {
        this->coeffs = (real_t*) malloc((other.order + 1) * sizeof(real_t));
        if (this->coeffs != nullptr) {
            this->order = other.order;
            memcpy(this->coeffs, other.coeffs, (other.order + 1) * sizeof(real_t));
        }
    }
}

// copy
void Polynomial::operator=(const Polynomial& other) {
    if (this != &other) {
        if (this->coeffs != nullptr) { free(this->coeffs); this->coeffs = nullptr; }
        if (other.coeffs != nullptr) {
            this->coeffs = (real_t*) malloc((order + 1) * sizeof(real_t));
            if (this->coeffs != nullptr) {
                this->order = other.order;
                memcpy(this->coeffs, other.coeffs, (order + 1) * sizeof(real_t));
            }
        }
    }
}

// move
Polynomial::Polynomial(Polynomial&& other):
order(other.order),
coeffs(other.coeffs)
{
    free(other.coeffs);
    other.coeffs = nullptr;
    other.order = 0;
}

// move
void Polynomial::operator=(Polynomial&& other) {
    if (this != &other) {
        if (this->coeffs != nullptr) { free(this->coeffs); this->coeffs = nullptr; }
        this->order = other.order;
        this->coeffs = other.coeffs;
        other.order = 0;
        other.coeffs = nullptr;
    }
}

// polyval
real_t Polynomial::operator()(const real_t x) const {
    if (this-> coeffs == nullptr) { return 0; }

    // Horner
    real_t accumulator = this->coeffs[this->order];
    for (size_t k = this->order; k > 0; /*--k*/) {
        accumulator = x * accumulator + this->coeffs[--k /*k-1*/] ;
    }

    // return
    return accumulator;
}

// access coeffs
real_t& Polynomial::operator[](const size_t index) {
    return this->coeffs[index];
}

// access coeffs
const real_t& Polynomial::operator[](const size_t index) const {
    return this->coeffs[index];
}

// access coeffs
real_t& Polynomial::a(const size_t index) {
    return this->coeffs[index];
}

// access coeffs
const real_t& Polynomial::a(const size_t index) const {
    return this->coeffs[index];
}

// convert to readable form
std::string Polynomial::to_string() const {

    // std::stringstream buff;
    // buff << "{";
    // for (size_t k = 0; k <= this->order; ++k) {
    //     if (k < this->order) { buff << this->a(k) << ", "; }
    //     else { buff << this->a(k); }
    // }
    // buff << "}";

    // return buff.str();

    std::stringstream buff;
    buff << "{";
    for (size_t k = 0; k <= this->order; ++k) {
        if (k < this->order) { buff << this->a(this->order - k) << ", "; }
        else { buff << this->a(this->order - k); }
    }
    buff << "}";

    return buff.str();

}

// Destructor
Polynomial::~Polynomial() {
    if (this->coeffs != nullptr) {
        free(this->coeffs);
    }
}

// Getter
size_t Polynomial::get_order() const {
    return this->order;
}

// operators
Polynomial Polynomial::operator+() const {
    // Copy
    Polynomial result(*this);
    result.shrink_to_fit();
    return result;
}

Polynomial Polynomial::operator-() const {
    // Copy
    Polynomial result(this->order);
    for (size_t k = 0; k < this->order; ++k) {
        result.a(k) = (-1) * this->a(k);
    }

    result.shrink_to_fit();
    return result;
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
    if (this->order >= other.order) {
        Polynomial result(this->order);

        for (size_t k = 0; k <= other.order; ++k) { 
            result.a(k) = this->a(k) + other.a(k);
        }

        for (size_t k = other.order + 1; k <= this->order; ++k) {
            result.a(k) = this->a(k);
        }

        result.shrink_to_fit();
        return result;
    }
    else {
        Polynomial result(this->order);

        for (size_t k = 0; k <= this->order; ++k) { 
            result.a(k) = this->a(k) + other.a(k);
        }

        for (size_t k = this->order + 1; k <= other.order; ++k) {
            result.a(k) = other.a(k);
        }

        result.shrink_to_fit();
        return result;
    }
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
    if (this->order >= other.order) {
        Polynomial result(this->order);

        for (size_t k = 0; k <= other.order; ++k) { 
            result.a(k) = this->a(k) - other.a(k);
        }

        for (size_t k = other.order + 1; k <= this->order; ++k) {
            result.a(k) = this->a(k);
        }

        result.shrink_to_fit();
        return result;
    }
    else {
        Polynomial result(this->order);

        for (size_t k = 0; k <= this->order; ++k) { 
            result.a(k) = this->a(k) - other.a(k);
        }

        for (size_t k = this->order + 1; k <= other.order; ++k) {
            result.a(k) = (-1) * other.a(k);
        }

        result.shrink_to_fit();
        return result;
    }
}

Polynomial Polynomial::operator*(const Polynomial& other) const {

    // From Matlab:
    // w = conv(u,v) returns the convolution of vectors u and v. 
    // If u and v are vectors of polynomial coefficients, 
    // convolving them is equivalent to multiplying the two polynomials.

    Polynomial result(this->order + other.order);
    for (size_t	k = 0; k <= result.order; k++) {

        const size_t start = (k <= this->order ? 0 : k - this->order);
        const size_t end = (k <= other.order ? k : other.order);
        result.a(k) = 0;

        for (size_t j = start; j <= end; ++j) {
            result.a(k) += other.a(j) * this->a(k-j);
        }
    }

    result.shrink_to_fit();
    return result;
}


void Polynomial::operator+=(const Polynomial& other) {
    Polynomial result = (*this) + other;
    (*this) = std::move(result);
}

void Polynomial::operator-=(const Polynomial& other) {
    Polynomial result = (*this) - other;
    (*this) = std::move(result);
}

void Polynomial::operator*=(const Polynomial& other) {
    Polynomial result = (*this) * other;
    (*this) = std::move(result);
}


Polynomial Polynomial::operator+(const real_t scalar) const {
    // Copy
    Polynomial result(*this);
    result.a(0) += scalar;

    result.shrink_to_fit();
    return result;
}

Polynomial Polynomial::operator-(const real_t scalar) const {
    // Copy
    Polynomial result(*this);
    result.a(0) -= scalar;

    result.shrink_to_fit();
    return result;
}

Polynomial Polynomial::operator*(const real_t scalar) const {
    Polynomial result(this->order);
    for (size_t k = 0; k <= this->order; ++k) {
        result.a(k) = this->a(k) * scalar ;
    }

    result.shrink_to_fit();
    return result;
}

Polynomial Polynomial::operator/(const real_t scalar) const {
    Polynomial result(this->order);
    for (size_t k = 0; k <= this->order; ++k) {
        result.a(k) = this->a(k) / scalar;
    }

    result.shrink_to_fit();
    return result;
}



void Polynomial::operator+=(const real_t scalar) {
    this->a(0) += scalar;
    this->shrink_to_fit();
}

void Polynomial::operator-=(const real_t scalar) {
    this->a(0) -= scalar;
    this->shrink_to_fit();
}

void Polynomial::operator*=(const real_t scalar) {
    this->shrink_to_fit();
    for (size_t k = 0; k <= this->order; ++k) {
        this->a(k) *= scalar;
    }
}

void Polynomial::operator/=(const real_t scalar) {
    this->shrink_to_fit();
    for (size_t k = 0; k <= this->order; ++k) {
        this->a(k) /= scalar;
    }
}


Polynomial operator+(const real_t scalar, const Polynomial& poly) {
    // Copy
    Polynomial result(poly);
    result.a(0) += scalar;

    result.shrink_to_fit();
    return result;
}

Polynomial operator-(const real_t scalar, const Polynomial& poly) {
    // Copy
    Polynomial result(poly);
    result.a(0) -= scalar;

    result.shrink_to_fit();
    return result;
}

Polynomial operator*(const real_t scalar, const Polynomial& poly) {
    Polynomial result(poly.get_order());
    for (size_t k = 0; k <= poly.get_order(); ++k) {
        result.a(k) = poly.a(k) * scalar ;
    }

    result.shrink_to_fit();
    return result;
}


// polydiv
std::pair<Polynomial, Polynomial> Polynomial::operator/(const Polynomial& divisor) const {
    this->shrink_to_fit();
    divisor.shrink_to_fit();

    if (this->order < divisor.order) {
        Polynomial quotient(0, NULL);
        Polynomial remainder(*this);
        return std::make_pair(quotient, remainder);
    }

    Polynomial quotient(this->order - divisor.order);
    Polynomial remainder(*this);

    size_t n = remainder.order;
    const size_t m = divisor.order;
    for (size_t k = quotient.order + 1; k > 0; /*--k*/) {
        
        // Do this at the start of the loop
        --k;
        
        // Divide
        quotient.a(k) = remainder.a(n) / divisor.a(m);

        // Subtract
        remainder.a(n) = 0;
        for (size_t j = 0; j < m; ++j) {
            remainder.a(j + k) -= quotient.a(k) * divisor.a(j);
        }

        // Reduce
        --n;
    }

    // quotient.shrink_to_fit();
    remainder.shrink_to_fit();
    return std::make_pair(quotient, remainder);
}

// polyfit
Polynomial Polynomial::fit(const size_t order, const std::vector<real_t>& x, const std::vector<real_t>& y) {

    const auto coeff_mat = Matrix::solve(Matrix::vandermonde(order, x), Matrix::column_vector(y));
    if (coeff_mat.get_elements() != nullptr) {
        return Polynomial(order, coeff_mat.get_elements());
    }
    else {
        return Polynomial();
    }
}
void Polynomial::fit(real_t* const p, const size_t order, const real_t* const x, const real_t* const y, const size_t size) {
    std::vector<real_t> x_vec(x, x+size);
    std::vector<real_t> y_vec(y, y+size);
    const auto poly = Polynomial::fit(order, x_vec, y_vec);
    memcpy(p, poly.coeffs, (order + 1) * sizeof(real_t));
}


void Polynomial::shrink_to_fit() const {

    // if (an != 0) => Nothing to do
    if ((this->order == 0) || (this->a(this->order) != 0)) { return; }

    // Don't check a0 => (k >= 1)
    for (size_t k = this->order-1; k >= 1; --k) {

        // find first nonzero element
        if ((this->a(k) != 0)) {

            // shrink array
            real_t* const new_array = (real_t*) realloc(this->coeffs, (k + 1) * sizeof(real_t));

            // check
            if (new_array != nullptr) { 
                this->order = k; 
                this->coeffs = new_array;
            }
            return;
        }
    }

    // shrink array to size 1 (order = 0)
    real_t* const new_array = (real_t*) realloc(this->coeffs, 1 * sizeof(real_t));

    // check
    if (new_array != nullptr) { 
        this->order = 0; 
        this->coeffs = new_array;
    }
}

void Polynomial::grow_to(const size_t new_order) const {
    if (new_order <= this->order) { return; }
    const size_t old_order = this->order;

    // grow array
    real_t* const new_array = (real_t*) realloc(this->coeffs, (new_order + 1) * sizeof(real_t));

    // check
    if (new_array != nullptr) { 

        // initilize new coeffs
        memset(&(new_array[old_order]), 0, (new_order - old_order) * sizeof(real_t));
        this->order = new_order; 
        this->coeffs = new_array;
    }
    return;
}


// ----- iterator -----

// Constructor
Polynomial::iterator::iterator(Polynomial& poly, const int32_t index):
poly(&poly),
index(index)
{}
Polynomial::iterator::iterator(Polynomial* const poly, const int32_t index):
poly(poly),
index(index)
{

}


// Comparison
bool Polynomial::iterator::operator==(const Polynomial::iterator& other) const {
    return this->index == other.index;
}
bool Polynomial::iterator::operator!=(const Polynomial::iterator& other) const {
    return this->index != other.index;
}
bool Polynomial::iterator::operator<(const Polynomial::iterator& other) const {
    return this->index < other.index;
}
bool Polynomial::iterator::operator<=(const Polynomial::iterator& other) const {
    return this->index <= other.index;
}
bool Polynomial::iterator::operator>(const Polynomial::iterator& other) const {
    return this->index > other.index;
}
bool Polynomial::iterator::operator>=(const Polynomial::iterator& other) const {
    return this->index >= other.index;
}

// Access
Polynomial::iterator::reference Polynomial::iterator::operator*() {
    return this->poly->a(index);
}
Polynomial::iterator::const_reference Polynomial::iterator::operator*() const {
    return this->poly->a(index);
}
Polynomial::iterator::reference Polynomial::iterator::operator->() {
    return this->poly->a(index);
}
Polynomial::iterator::const_reference Polynomial::iterator::operator->() const {
    return this->poly->a(index);
}
Polynomial::iterator::reference Polynomial::iterator::operator[](const int32_t offset) {
    return this->poly->a(offset + index);
}
Polynomial::iterator::const_reference Polynomial::iterator::operator[](const int32_t offset)const {
    return this->poly->a(offset + index);
}

// Prefix increment
Polynomial::iterator& Polynomial::iterator::operator++() {
    this->index += 1;
    return *this;
}
Polynomial::iterator& Polynomial::iterator::operator--() {
    this->index -= 1;
    return *this;
}

// Postfix increment
Polynomial::iterator Polynomial::iterator::operator++(int) {
    iterator other(*this);
    this->index += 1;
    return other;
}
Polynomial::iterator Polynomial::iterator::operator--(int) {
    iterator other(*this);
    this->index -= 1;
    return other;
}

// Arithmetic
Polynomial::iterator Polynomial::iterator::operator+(const int32_t n) const {
    return iterator(this->poly, this->index + n);
}
Polynomial::iterator Polynomial::iterator::operator-(const int32_t n) const {
    return iterator(this->poly, this->index - n);
}
Polynomial::iterator& Polynomial::iterator::operator+=(const int32_t n) {
    this->index += n;
    return *this;
}
Polynomial::iterator& Polynomial::iterator::operator-=(const int32_t n) {
    this->index -= n;
    return *this;
}
int32_t Polynomial::iterator::operator-(const Polynomial::iterator& other) const {
    return (this->index - other.index);
}





// ----- const_iterator -----

// Constructor
Polynomial::const_iterator::const_iterator(const Polynomial& poly, const int32_t index):
poly(&poly),
index(index)
{}
Polynomial::const_iterator::const_iterator(const Polynomial* const poly, const int32_t index):
poly(poly),
index(index)
{

}


// Comparison
bool Polynomial::const_iterator::operator==(const Polynomial::const_iterator& other) const {
    return this->index == other.index;
}
bool Polynomial::const_iterator::operator!=(const Polynomial::const_iterator& other) const {
    return this->index != other.index;
}
bool Polynomial::const_iterator::operator<(const Polynomial::const_iterator& other) const {
    return this->index < other.index;
}
bool Polynomial::const_iterator::operator<=(const Polynomial::const_iterator& other) const {
    return this->index <= other.index;
}
bool Polynomial::const_iterator::operator>(const Polynomial::const_iterator& other) const {
    return this->index > other.index;
}
bool Polynomial::const_iterator::operator>=(const Polynomial::const_iterator& other) const {
    return this->index >= other.index;
}

// Access
Polynomial::const_iterator::const_reference Polynomial::const_iterator::operator*() const {
    return this->poly->a(index);
}
Polynomial::const_iterator::const_reference Polynomial::const_iterator::operator->() const {
    return this->poly->a(index);
}
Polynomial::const_iterator::const_reference Polynomial::const_iterator::operator[](const int32_t offset)const {
    return this->poly->a(offset + index);
}

// Prefix increment
Polynomial::const_iterator& Polynomial::const_iterator::operator++() {
    this->index += 1;
    return *this;
}
Polynomial::const_iterator& Polynomial::const_iterator::operator--() {
    this->index -= 1;
    return *this;
}

// Postfix increment
Polynomial::const_iterator Polynomial::const_iterator::operator++(int) {
    const_iterator other(*this);
    this->index += 1;
    return other;
}
Polynomial::const_iterator Polynomial::const_iterator::operator--(int) {
    const_iterator other(*this);
    this->index -= 1;
    return other;
}

// Arithmetic
Polynomial::const_iterator Polynomial::const_iterator::operator+(const int32_t n) const {
    return const_iterator(this->poly, this->index + n);
}
Polynomial::const_iterator Polynomial::const_iterator::operator-(const int32_t n) const {
    return const_iterator(this->poly, this->index - n);
}
Polynomial::const_iterator& Polynomial::const_iterator::operator+=(const int32_t n) {
    this->index += n;
    return *this;
}
Polynomial::const_iterator& Polynomial::const_iterator::operator-=(const int32_t n) {
    this->index -= n;
    return *this;
}
int32_t Polynomial::const_iterator::operator-(const Polynomial::const_iterator& other) const {
    return (this->index - other.index);
}



// ----- reverse_iterator -----

// Constructor
Polynomial::reverse_iterator::reverse_iterator(Polynomial& poly, const int32_t index):
poly(&poly),
index(index)
{}
Polynomial::reverse_iterator::reverse_iterator(Polynomial* const poly, const int32_t index):
poly(poly),
index(index)
{

}


// Comparison
bool Polynomial::reverse_iterator::operator==(const Polynomial::reverse_iterator& other) const {
    return this->index == other.index;
}
bool Polynomial::reverse_iterator::operator!=(const Polynomial::reverse_iterator& other) const {
    return this->index != other.index;
}
bool Polynomial::reverse_iterator::operator<(const Polynomial::reverse_iterator& other) const {
    return this->index < other.index;
}
bool Polynomial::reverse_iterator::operator<=(const Polynomial::reverse_iterator& other) const {
    return this->index <= other.index;
}
bool Polynomial::reverse_iterator::operator>(const Polynomial::reverse_iterator& other) const {
    return this->index > other.index;
}
bool Polynomial::reverse_iterator::operator>=(const Polynomial::reverse_iterator& other) const {
    return this->index >= other.index;
}

// Access
Polynomial::reverse_iterator::reference Polynomial::reverse_iterator::operator*() {
    return this->poly->a(this->poly->get_order() - index);
}
const Polynomial::reverse_iterator::reference Polynomial::reverse_iterator::operator*() const {
    return this->poly->a(this->poly->get_order() - index);
}
Polynomial::reverse_iterator::reference Polynomial::reverse_iterator::operator->() {
    return this->poly->a(this->poly->get_order() - index);
}
const Polynomial::reverse_iterator::reference Polynomial::reverse_iterator::operator->() const {
    return this->poly->a(this->poly->get_order() - index);
}
Polynomial::reverse_iterator::reference Polynomial::reverse_iterator::operator[](const int32_t offset) {
    return this->poly->a(this->poly->get_order() - (offset + index));
}
const Polynomial::reverse_iterator::reference Polynomial::reverse_iterator::operator[](const int32_t offset)const {
    return this->poly->a(this->poly->get_order() - (offset + index));
}

// Prefix increment
Polynomial::reverse_iterator& Polynomial::reverse_iterator::operator++() {
    this->index += 1;
    return *this;
}
Polynomial::reverse_iterator& Polynomial::reverse_iterator::operator--() {
    this->index -= 1;
    return *this;
}

// Postfix increment
Polynomial::reverse_iterator Polynomial::reverse_iterator::operator++(int) {
    reverse_iterator other(*this);
    this->index += 1;
    return other;
}
Polynomial::reverse_iterator Polynomial::reverse_iterator::operator--(int) {
    reverse_iterator other(*this);
    this->index -= 1;
    return other;
}

// Arithmetic
Polynomial::reverse_iterator Polynomial::reverse_iterator::operator+(const int32_t n) const {
    return reverse_iterator(this->poly, this->index + n);
}
Polynomial::reverse_iterator Polynomial::reverse_iterator::operator-(const int32_t n) const {
    return reverse_iterator(this->poly, this->index - n);
}
Polynomial::reverse_iterator& Polynomial::reverse_iterator::operator+=(const int32_t n) {
    this->index += n;
    return *this;
}
Polynomial::reverse_iterator& Polynomial::reverse_iterator::operator-=(const int32_t n) {
    this->index -= n;
    return *this;
}
int32_t Polynomial::reverse_iterator::operator-(const Polynomial::reverse_iterator& other) const {
    return (this->index - other.index);
}




// ----- const_reverse_iterator -----

// Constructor
Polynomial::const_reverse_iterator::const_reverse_iterator(const Polynomial& poly, const int32_t index):
poly(&poly),
index(index)
{}
Polynomial::const_reverse_iterator::const_reverse_iterator(const Polynomial* const poly, const int32_t index):
poly(poly),
index(index)
{

}


// Comparison
bool Polynomial::const_reverse_iterator::operator==(const Polynomial::const_reverse_iterator& other) const {
    return this->index == other.index;
}
bool Polynomial::const_reverse_iterator::operator!=(const Polynomial::const_reverse_iterator& other) const {
    return this->index != other.index;
}
bool Polynomial::const_reverse_iterator::operator<(const Polynomial::const_reverse_iterator& other) const {
    return this->index < other.index;
}
bool Polynomial::const_reverse_iterator::operator<=(const Polynomial::const_reverse_iterator& other) const {
    return this->index <= other.index;
}
bool Polynomial::const_reverse_iterator::operator>(const Polynomial::const_reverse_iterator& other) const {
    return this->index > other.index;
}
bool Polynomial::const_reverse_iterator::operator>=(const Polynomial::const_reverse_iterator& other) const {
    return this->index >= other.index;
}

// Access
Polynomial::const_reverse_iterator::const_reference Polynomial::const_reverse_iterator::operator*() const {
    return this->poly->a(this->poly->get_order() - index);
}
Polynomial::const_reverse_iterator::const_reference Polynomial::const_reverse_iterator::operator->() const {
    return this->poly->a(this->poly->get_order() - index);
}
Polynomial::const_reverse_iterator::const_reference Polynomial::const_reverse_iterator::operator[](const int32_t offset)const {
    return this->poly->a(this->poly->get_order() - (offset + index));
}

// Prefix increment
Polynomial::const_reverse_iterator& Polynomial::const_reverse_iterator::operator++() {
    this->index += 1;
    return *this;
}
Polynomial::const_reverse_iterator& Polynomial::const_reverse_iterator::operator--() {
    this->index -= 1;
    return *this;
}

// Postfix increment
Polynomial::const_reverse_iterator Polynomial::const_reverse_iterator::operator++(int) {
    const_reverse_iterator other(*this);
    this->index += 1;
    return other;
}
Polynomial::const_reverse_iterator Polynomial::const_reverse_iterator::operator--(int) {
    const_reverse_iterator other(*this);
    this->index -= 1;
    return other;
}

// Arithmetic
Polynomial::const_reverse_iterator Polynomial::const_reverse_iterator::operator+(const int32_t n) const {
    return const_reverse_iterator(this->poly, this->index + n);
}
Polynomial::const_reverse_iterator Polynomial::const_reverse_iterator::operator-(const int32_t n) const {
    return const_reverse_iterator(this->poly, this->index - n);
}
Polynomial::const_reverse_iterator& Polynomial::const_reverse_iterator::operator+=(const int32_t n) {
    this->index += n;
    return *this;
}
Polynomial::const_reverse_iterator& Polynomial::const_reverse_iterator::operator-=(const int32_t n) {
    this->index -= n;
    return *this;
}
int32_t Polynomial::const_reverse_iterator::operator-(const Polynomial::const_reverse_iterator& other) const {
    return (this->index - other.index);
}




Polynomial::iterator Polynomial::begin() {
    return iterator(this, 0);
}
Polynomial::const_iterator Polynomial::begin() const {
    return const_iterator(this, 0);
}
Polynomial::iterator Polynomial::end() {
    return iterator(this, this->order+1);
}
Polynomial::const_iterator Polynomial::end() const {
    return const_iterator(this, this->order+1);
}
Polynomial::reverse_iterator Polynomial::rbegin() {
    return reverse_iterator(this, 0);
}
Polynomial::const_reverse_iterator Polynomial::rbegin() const {
    return const_reverse_iterator(this, 0);
}
Polynomial::reverse_iterator Polynomial::rend() {
    return reverse_iterator(this, this->order+1);
}
Polynomial::const_reverse_iterator Polynomial::rend() const {
    return const_reverse_iterator(this, this->order+1);
}
Polynomial::const_iterator Polynomial::cbegin() const {
    return const_iterator(this, 0);
}
Polynomial::const_iterator Polynomial::cend() const {
    return const_iterator(this, this->order+1);
}
Polynomial::const_reverse_iterator Polynomial::crbegin() const {
    return const_reverse_iterator(this, 0);
}
Polynomial::const_reverse_iterator Polynomial::crend() const {
    return const_reverse_iterator(this, this->order+1);
}

