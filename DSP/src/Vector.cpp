#include <cmath> // absf, acosf, fsqrt
#include "DSP/Math/Vector.hpp"

// #include <utility> // move


using namespace DSP::Math;


// Constructor
Vector::Vector():
Matrix()
{}
Vector::Vector(const size_t size):
Matrix(size, 1)
{}
Vector::Vector(const size_t size, const real_t* vec):
Matrix(size, 1, vec)
{}
Vector::Vector(const std::vector<real_t> vec):
Matrix(Matrix::column_vector(vec))
{}

// Initilize
Vector::Vector(const std::initializer_list<real_t>& il):
Matrix(Matrix(il).transpose())
{}
void Vector::operator=(const std::initializer_list<real_t>& il) {
    this->Matrix::operator=(Matrix(il).transpose());
}

// Copy
Vector::Vector(const Matrix& column_mat):
Matrix()
{
    if (column_mat.get_num_columns() == 1) {
        this->Matrix::operator=(column_mat);
    }
}
Vector::Vector(const Vector& other):
Matrix(other)
{}
void Vector::operator=(const Vector& other) {
    this->Matrix::operator=(other);
}

// Move
Vector::Vector(Matrix&& column_mat):
Matrix()
{
    if (column_mat.get_num_columns() == 1) {
        this->Matrix::operator=(column_mat);
    }
}
Vector::Vector(Vector&& other):
Matrix(std::move(other))
{}
void Vector::operator=(Vector&& other) {
    this->Matrix::operator=(std::move(other));
}

// Destructor
Vector::~Vector() {
    // Noting to do
}

// to Matrix
Matrix Vector::to_Matrix() const {
    return Matrix(this->rows, this->columns, this->elements);
}


// to real_t
Vector::operator real_t() const {
    return this->length();
}
real_t Vector::abs() const {
    return Vector::abs(*this);
}
real_t Vector::norm_squared() const {
    return Vector::norm_squared(*this);
}
real_t Vector::norm() const {
    return Vector::norm(*this);
}
real_t Vector::length_squared() const {
    return Vector::length_squared(*this);
}
real_t Vector::length() const {
    return Vector::length(*this);
}
Vector Vector::normalized() const {
    return Vector::normalized(*this);
}
real_t Vector::cosphi(const Vector& other) const {
    return Vector::cosphi(*this, other);
}
real_t Vector::angle(const Vector& other) const {
    return Vector::angle(*this, other);
}
real_t Vector::project(const Vector& other) const {
    return Vector::project(*this, other);
}
Vector Vector::projection(const Vector& other) const {
    return Vector::projection(*this, other);
}




// element access 1-based
real_t& Vector::operator()(const size_t num_element) {
    return this->elements[num_element - 1];
}
const real_t& Vector::operator()(const size_t num_element) const {
    return this->elements[num_element - 1];
}

// Element access 0-based
real_t& Vector::v(const size_t index) {
    return this->elements[index];
}
const real_t& Vector::v(const size_t index) const {
    return this->elements[index];
}

// Element access 1-based
real_t& Vector::at(const size_t num_element) {
    return this->elements[num_element - 1];
}
const real_t& Vector::at(const size_t num_element) const {
    return this->elements[num_element - 1];
}



// Arithmetic
Vector Vector::operator+() const {
    return *this;
}
Vector Vector::operator-() const {
    Vector other(this->size());
    for (size_t k = 0; k < this->size(); ++k) {
        other.elements[k] = -this->elements[k];
    }
    return other;
}

// Arithmetic
Vector Vector::operator+(const Vector& other) const {
    return this->Matrix::operator+(other);
}
Vector Vector::operator-(const Vector& other) const {
    return this->Matrix::operator-(other);
}
real_t Vector::operator*(const Vector& other) const {
    return this->dot(other);
}

// Arithmetic
void Vector::operator+=(const Vector& other) {
    this->Matrix::operator+=(other);
}
void Vector::operator-=(const Vector& other) {
    this->Matrix::operator-=(other);
}

// Arithmetic
Vector Vector::operator*(const real_t scalar) const {
    return this->Matrix::operator*(scalar);
}
Vector Vector::operator/(const real_t scalar) const {
    return this->Matrix::operator/(scalar);
}
void Vector::operator*=(const real_t scalar) {
    this->Matrix::operator*=(scalar);
}
void Vector::operator/=(const real_t scalar) {
    this->Matrix::operator/=(scalar);
}



real_t Vector::dot(const Vector& other) const {
    return Vector::dot(*this, other);
}
Vector Vector::cross(const Vector& other) const {
    return Vector::cross(*this, other);
}
Vector Vector::conv(const Vector& other) const {
    return Vector::conv(*this, other);
}
std::pair<Vector, Vector> Vector::deconv(const Vector& other) const {
    return Vector::deconv(*this, other);
}


real_t Vector::dot(const Vector& u, const Vector& v) {
    if (u.size() != v.size()) { return 0; }
    if (u.size() == 0) { return 0; }

    real_t sum = u.elements[0] * v.elements[0];
    for (size_t k = 0; k < u.size(); ++k) {
        sum += u.elements[k] * v.elements[k];
    }

    return sum;
}
Vector Vector::cross(const Vector& u, const Vector& v) {
    if (u.size() != 3 && v.size() != 3) { return Vector(); }
    
    const Vector w = {
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0]
    };
    return w;
}
Vector Vector::conv(const Vector& u, const Vector& v) {
    //Check
    if (u.size() == 0 || v.size() == 0) { return Vector(); }
    const size_t conv_size = (u.size() + v.size() - 1);

    // From Matlab:
    // w = conv(u,v) returns the convolution of vectors u and v. 
    // If u and v are vectors of polynomial coefficients, 
    // convolving them is equivalent to multiplying the two polynomials.
    size_t start, end;
    Vector w(conv_size);
    for (size_t	k = 0; k < conv_size; k++) {

        // Determine start and end
        start = (k < u.size() ? 0 : k - (u.size() - 1));
        end = (k < v.size() ? k : (v.size() - 1));

        // Init
        w[k] = 0;

        // Calculate k-th element
        for (size_t j = start; j <= end; ++j) {
            w[k] += v[j] * u[k-j];
        }
    }

    return w;
}
std::pair<Vector, Vector> Vector::deconv(const Vector& u, const Vector& v) {

    // Check
    if (v[0] == 0) { return std::make_pair(Vector(), Vector()); }
    if (u.size() < v.size()) { return std::make_pair(Vector(), Vector());; }        

    // calculate sizes
    const size_t q_size = u.size() - v.size() + 1;
    // const size_t r_size = v.size() - 1;

    // we copy 'u' into 'r' at the start
    Vector r = u;
    Vector q(q_size);

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



real_t Vector::abs(const Vector& vec) {
    Vector vec2(vec.size());
    for (size_t k = 0; k < vec.size(); ++k) {
        vec2.elements[k] = (vec.elements[k] >= 0 ? vec.elements[k] : -(vec.elements[k]));
    }
    return vec2;
}
real_t Vector::norm_squared(const Vector& vec) {
    if (vec.size() == 0) { return 0; }
    real_t sum = vec.elements[0] * vec.elements[0];
    for (size_t k = 0; k < vec.size(); ++k) {
        sum += vec.elements[k] * vec.elements[k];
    }
    return sum;
}
real_t Vector::norm(const Vector& vec) {
    return sqrtf(Vector::norm_squared(vec));
}
real_t Vector::length_squared(const Vector& vec) {
    return Vector::norm_squared(vec);
}
real_t Vector::length(const Vector& vec) {
    return Vector::norm(vec);
}
Vector Vector::normalized(const Vector& vec) {
    const real_t vec_length = vec.length();
    Vector vec2(vec.size());
    for (size_t k = 0; k < vec.size(); ++k) {
        vec2.elements[k] = vec.elements[k] / vec_length;
    }
    return vec2;
}
real_t Vector::cosphi(const Vector& a, const Vector& b) {
    return (a * b) / (a.length() * b.length());
}
real_t Vector::angle(const Vector& a, const Vector& b) {
    return acosf(Vector::cosphi(a, b));
}
real_t Vector::project(const Vector& a, const Vector& b) {
    return (a * b) / b.length();
}
Vector Vector::projection(const Vector& a, const Vector& b) {
    // return Vector::normalized(b) * Vector::project(a, b);
    return ((a * b) / b.length_squared()) * b;
}


Vector operator*(const real_t scalar, const Vector& v) {
    return v * scalar;
}
Vector operator*(const Matrix& M, const Vector& v) {
    // // return Vector(M.operator*(v));
    // return Vector(M * v.to_Matrix());
    Vector w(M.get_num_rows(), nullptr);
    for (size_t k = 0; k < w.size(); ++k) {
        for (size_t j = 0; j < v.size(); ++j) {
            w.v(k) += M.m(k, j) * v.v(j);
        }
    }
    return w;
}









// ----- iterator -----

// Constructor
Vector::iterator::iterator(Vector& vec, const int32_t index):
vec(&vec),
index(index)
{}
Vector::iterator::iterator(Vector* const vec, const int32_t index):
vec(vec),
index(index)
{

}


// Comparison
bool Vector::iterator::operator==(const Vector::iterator& other) const {
    return this->index == other.index;
}
bool Vector::iterator::operator!=(const Vector::iterator& other) const {
    return this->index != other.index;
}
bool Vector::iterator::operator<(const Vector::iterator& other) const {
    return this->index < other.index;
}
bool Vector::iterator::operator<=(const Vector::iterator& other) const {
    return this->index <= other.index;
}
bool Vector::iterator::operator>(const Vector::iterator& other) const {
    return this->index > other.index;
}
bool Vector::iterator::operator>=(const Vector::iterator& other) const {
    return this->index >= other.index;
}

// Access
Vector::iterator::reference Vector::iterator::operator*() {
    return this->vec->v(index);
}
Vector::iterator::const_reference Vector::iterator::operator*() const {
    return this->vec->v(index);
}
Vector::iterator::reference Vector::iterator::operator->() {
    return this->vec->v(index);
}
Vector::iterator::const_reference Vector::iterator::operator->() const {
    return this->vec->v(index);
}
Vector::iterator::reference Vector::iterator::operator[](const int32_t offset) {
    return this->vec->v(offset + index);
}
Vector::iterator::const_reference Vector::iterator::operator[](const int32_t offset)const {
    return this->vec->v(offset + index);
}

// Prefix increment
Vector::iterator& Vector::iterator::operator++() {
    this->index += 1;
    return *this;
}
Vector::iterator& Vector::iterator::operator--() {
    this->index -= 1;
    return *this;
}

// Postfix increment
Vector::iterator Vector::iterator::operator++(int) {
    iterator other(*this);
    this->index += 1;
    return other;
}
Vector::iterator Vector::iterator::operator--(int) {
    iterator other(*this);
    this->index -= 1;
    return other;
}

// Arithmetic
Vector::iterator Vector::iterator::operator+(const int32_t n) const {
    return iterator(this->vec, this->index + n);
}
Vector::iterator Vector::iterator::operator-(const int32_t n) const {
    return iterator(this->vec, this->index - n);
}
Vector::iterator& Vector::iterator::operator+=(const int32_t n) {
    this->index += n;
    return *this;
}
Vector::iterator& Vector::iterator::operator-=(const int32_t n) {
    this->index -= n;
    return *this;
}
int32_t Vector::iterator::operator-(const Vector::iterator& other) const {
    return (this->index - other.index);
}





// ----- const_iterator -----

// Constructor
Vector::const_iterator::const_iterator(const Vector& vec, const int32_t index):
vec(&vec),
index(index)
{}
Vector::const_iterator::const_iterator(const Vector* const vec, const int32_t index):
vec(vec),
index(index)
{

}


// Comparison
bool Vector::const_iterator::operator==(const Vector::const_iterator& other) const {
    return this->index == other.index;
}
bool Vector::const_iterator::operator!=(const Vector::const_iterator& other) const {
    return this->index != other.index;
}
bool Vector::const_iterator::operator<(const Vector::const_iterator& other) const {
    return this->index < other.index;
}
bool Vector::const_iterator::operator<=(const Vector::const_iterator& other) const {
    return this->index <= other.index;
}
bool Vector::const_iterator::operator>(const Vector::const_iterator& other) const {
    return this->index > other.index;
}
bool Vector::const_iterator::operator>=(const Vector::const_iterator& other) const {
    return this->index >= other.index;
}

// Access
Vector::const_iterator::const_reference Vector::const_iterator::operator*() const {
    return this->vec->v(index);
}
Vector::const_iterator::const_reference Vector::const_iterator::operator->() const {
    return this->vec->v(index);
}
Vector::const_iterator::const_reference Vector::const_iterator::operator[](const int32_t offset)const {
    return this->vec->v(offset + index);
}

// Prefix increment
Vector::const_iterator& Vector::const_iterator::operator++() {
    this->index += 1;
    return *this;
}
Vector::const_iterator& Vector::const_iterator::operator--() {
    this->index -= 1;
    return *this;
}

// Postfix increment
Vector::const_iterator Vector::const_iterator::operator++(int) {
    const_iterator other(*this);
    this->index += 1;
    return other;
}
Vector::const_iterator Vector::const_iterator::operator--(int) {
    const_iterator other(*this);
    this->index -= 1;
    return other;
}

// Arithmetic
Vector::const_iterator Vector::const_iterator::operator+(const int32_t n) const {
    return const_iterator(this->vec, this->index + n);
}
Vector::const_iterator Vector::const_iterator::operator-(const int32_t n) const {
    return const_iterator(this->vec, this->index - n);
}
Vector::const_iterator& Vector::const_iterator::operator+=(const int32_t n) {
    this->index += n;
    return *this;
}
Vector::const_iterator& Vector::const_iterator::operator-=(const int32_t n) {
    this->index -= n;
    return *this;
}
int32_t Vector::const_iterator::operator-(const Vector::const_iterator& other) const {
    return (this->index - other.index);
}



// ----- reverse_iterator -----

// Constructor
Vector::reverse_iterator::reverse_iterator(Vector& vec, const int32_t index):
vec(&vec),
index(index)
{}
Vector::reverse_iterator::reverse_iterator(Vector* const vec, const int32_t index):
vec(vec),
index(index)
{

}


// Comparison
bool Vector::reverse_iterator::operator==(const Vector::reverse_iterator& other) const {
    return this->index == other.index;
}
bool Vector::reverse_iterator::operator!=(const Vector::reverse_iterator& other) const {
    return this->index != other.index;
}
bool Vector::reverse_iterator::operator<(const Vector::reverse_iterator& other) const {
    return this->index < other.index;
}
bool Vector::reverse_iterator::operator<=(const Vector::reverse_iterator& other) const {
    return this->index <= other.index;
}
bool Vector::reverse_iterator::operator>(const Vector::reverse_iterator& other) const {
    return this->index > other.index;
}
bool Vector::reverse_iterator::operator>=(const Vector::reverse_iterator& other) const {
    return this->index >= other.index;
}

// Access
Vector::reverse_iterator::reference Vector::reverse_iterator::operator*() {
    return this->vec->v((this->vec->size() - 1) - index);
}
const Vector::reverse_iterator::reference Vector::reverse_iterator::operator*() const {
    return this->vec->v((this->vec->size() - 1) - index);
}
Vector::reverse_iterator::reference Vector::reverse_iterator::operator->() {
    return this->vec->v((this->vec->size() - 1) - index);
}
const Vector::reverse_iterator::reference Vector::reverse_iterator::operator->() const {
    return this->vec->v((this->vec->size() - 1) - index);
}
Vector::reverse_iterator::reference Vector::reverse_iterator::operator[](const int32_t offset) {
    return this->vec->v((this->vec->size() - 1) - (offset + index));
}
const Vector::reverse_iterator::reference Vector::reverse_iterator::operator[](const int32_t offset)const {
    return this->vec->v((this->vec->size() - 1) - (offset + index));
}

// Prefix increment
Vector::reverse_iterator& Vector::reverse_iterator::operator++() {
    this->index += 1;
    return *this;
}
Vector::reverse_iterator& Vector::reverse_iterator::operator--() {
    this->index -= 1;
    return *this;
}

// Postfix increment
Vector::reverse_iterator Vector::reverse_iterator::operator++(int) {
    reverse_iterator other(*this);
    this->index += 1;
    return other;
}
Vector::reverse_iterator Vector::reverse_iterator::operator--(int) {
    reverse_iterator other(*this);
    this->index -= 1;
    return other;
}

// Arithmetic
Vector::reverse_iterator Vector::reverse_iterator::operator+(const int32_t n) const {
    return reverse_iterator(this->vec, this->index + n);
}
Vector::reverse_iterator Vector::reverse_iterator::operator-(const int32_t n) const {
    return reverse_iterator(this->vec, this->index - n);
}
Vector::reverse_iterator& Vector::reverse_iterator::operator+=(const int32_t n) {
    this->index += n;
    return *this;
}
Vector::reverse_iterator& Vector::reverse_iterator::operator-=(const int32_t n) {
    this->index -= n;
    return *this;
}
int32_t Vector::reverse_iterator::operator-(const Vector::reverse_iterator& other) const {
    return (this->index - other.index);
}




// ----- const_reverse_iterator -----

// Constructor
Vector::const_reverse_iterator::const_reverse_iterator(const Vector& vec, const int32_t index):
vec(&vec),
index(index)
{}
Vector::const_reverse_iterator::const_reverse_iterator(const Vector* const vec, const int32_t index):
vec(vec),
index(index)
{

}


// Comparison
bool Vector::const_reverse_iterator::operator==(const Vector::const_reverse_iterator& other) const {
    return this->index == other.index;
}
bool Vector::const_reverse_iterator::operator!=(const Vector::const_reverse_iterator& other) const {
    return this->index != other.index;
}
bool Vector::const_reverse_iterator::operator<(const Vector::const_reverse_iterator& other) const {
    return this->index < other.index;
}
bool Vector::const_reverse_iterator::operator<=(const Vector::const_reverse_iterator& other) const {
    return this->index <= other.index;
}
bool Vector::const_reverse_iterator::operator>(const Vector::const_reverse_iterator& other) const {
    return this->index > other.index;
}
bool Vector::const_reverse_iterator::operator>=(const Vector::const_reverse_iterator& other) const {
    return this->index >= other.index;
}

// Access
Vector::const_reverse_iterator::const_reference Vector::const_reverse_iterator::operator*() const {
    return this->vec->v((this->vec->size() - 1) - index);
}
Vector::const_reverse_iterator::const_reference Vector::const_reverse_iterator::operator->() const {
    return this->vec->v((this->vec->size() - 1) - index);
}
Vector::const_reverse_iterator::const_reference Vector::const_reverse_iterator::operator[](const int32_t offset)const {
    return this->vec->v((this->vec->size() - 1) - (offset + index));
}

// Prefix increment
Vector::const_reverse_iterator& Vector::const_reverse_iterator::operator++() {
    this->index += 1;
    return *this;
}
Vector::const_reverse_iterator& Vector::const_reverse_iterator::operator--() {
    this->index -= 1;
    return *this;
}

// Postfix increment
Vector::const_reverse_iterator Vector::const_reverse_iterator::operator++(int) {
    const_reverse_iterator other(*this);
    this->index += 1;
    return other;
}
Vector::const_reverse_iterator Vector::const_reverse_iterator::operator--(int) {
    const_reverse_iterator other(*this);
    this->index -= 1;
    return other;
}

// Arithmetic
Vector::const_reverse_iterator Vector::const_reverse_iterator::operator+(const int32_t n) const {
    return const_reverse_iterator(this->vec, this->index + n);
}
Vector::const_reverse_iterator Vector::const_reverse_iterator::operator-(const int32_t n) const {
    return const_reverse_iterator(this->vec, this->index - n);
}
Vector::const_reverse_iterator& Vector::const_reverse_iterator::operator+=(const int32_t n) {
    this->index += n;
    return *this;
}
Vector::const_reverse_iterator& Vector::const_reverse_iterator::operator-=(const int32_t n) {
    this->index -= n;
    return *this;
}
int32_t Vector::const_reverse_iterator::operator-(const Vector::const_reverse_iterator& other) const {
    return (this->index - other.index);
}




Vector::iterator Vector::begin() {
    return iterator(this, 0);
}
Vector::const_iterator Vector::begin() const {
    return const_iterator(this, 0);
}
Vector::iterator Vector::end() {
    return iterator(this, this->size());
}
Vector::const_iterator Vector::end() const {
    return const_iterator(this, this->size());
}
Vector::reverse_iterator Vector::rbegin() {
    return reverse_iterator(this, 0);
}
Vector::const_reverse_iterator Vector::rbegin() const {
    return const_reverse_iterator(this, 0);
}
Vector::reverse_iterator Vector::rend() {
    return reverse_iterator(this, this->size());
}
Vector::const_reverse_iterator Vector::rend() const {
    return const_reverse_iterator(this, this->size());
}
Vector::const_iterator Vector::cbegin() const {
    return const_iterator(this, 0);
}
Vector::const_iterator Vector::cend() const {
    return const_iterator(this, this->size());
}
Vector::const_reverse_iterator Vector::crbegin() const {
    return const_reverse_iterator(this, 0);
}
Vector::const_reverse_iterator Vector::crend() const {
    return const_reverse_iterator(this, this->size());
}

