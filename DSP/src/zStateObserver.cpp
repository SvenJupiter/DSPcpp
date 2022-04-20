#include "DSP/Discrete/zStateObserver.hpp"

using namespace DSP;
using namespace DSP::Discrete;

// Constructor
zStateObserver::zStateObserver():
nx(), nu(), ny(),
A(), B(), C(), D(), L(),
x()
{}

// Constructor
zStateObserver::zStateObserver(const size_t num_states, const size_t num_inputs, const size_t num_outputs):
nx(num_states), nu(num_inputs), ny(num_outputs),
A(num_states, num_states), B(num_states, num_inputs), C(num_outputs, num_states), D(num_outputs, num_inputs), L(num_states, num_outputs),
x(num_states)
{}

// Constructor
zStateObserver::zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& L):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(A.get_num_rows()),
A(A), B(B), 
C(Math::Matrix::eye(A.get_num_rows())),
D(Math::Matrix::eye(A.get_num_rows(), B.get_num_columns())),
L(L), 
x(A.get_num_rows(), nullptr)
{}

// Constructor
zStateObserver::zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& L):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(C.get_num_rows()),
A(A), B(B), C(C),
D(Math::Matrix::eye(C.get_num_rows(), B.get_num_columns())),
L(L),
x(A.get_num_rows(), nullptr)
{}

// Constructor
zStateObserver::zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& D, const Math::Matrix& L):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(C.get_num_rows()),
A(A), B(B), C(C), D(D), L(L),
x(A.get_num_rows(), nullptr)
{}

// Constructor
zStateObserver::zStateObserver(const Math::Polynomial& num, const Math::Polynomial& den, const Math::Matrix& L):
zStateObserver(zStateSpace(num, den), L)
{}

// Constructor
zStateObserver::zStateObserver(const zStateSpace& zss, const Math::Matrix& L):
nx(zss.get_nx()), nu(zss.get_nu()), ny(zss.get_ny()),
A(zss.get_A()),
B(zss.get_B()),
C(zss.get_C()),
D(zss.get_D()),
L(L),
x(zss.get_nx(), nullptr)
{}

// Constructor
zStateObserver::zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& L, const Math::Vector& x0):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(A.get_num_rows()),
A(A), B(B), 
C(Math::Matrix::eye(A.get_num_rows())),
D(Math::Matrix::eye(A.get_num_rows(), B.get_num_columns())),
L(L), 
x(x0)
{}

// Constructor
zStateObserver::zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& L, const Math::Vector& x0):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(C.get_num_rows()),
A(A), B(B), C(C),
D(Math::Matrix::eye(C.get_num_rows(), B.get_num_columns())),
L(L),
x(x0)
{}

// Constructor
zStateObserver::zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& D, const Math::Matrix& L, const Math::Vector& x0):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(C.get_num_rows()),
A(A), B(B), C(C), D(D), L(L),
x(x0)
{}

// Constructor
zStateObserver::zStateObserver(const Math::Polynomial& num, const Math::Polynomial& den, const Math::Matrix& L, const Math::Vector& x0):
zStateObserver(zStateSpace(num, den), L, x0)
{}

// Constructor
zStateObserver::zStateObserver(const zStateSpace& zss, const Math::Matrix& L, const Math::Vector& x0):
nx(zss.get_nx()), nu(zss.get_nu()), ny(zss.get_ny()),
A(zss.get_A()),
B(zss.get_B()),
C(zss.get_C()),
D(zss.get_D()),
L(L),
x(x0)
{}




// Copy
zStateObserver::zStateObserver(const zStateObserver& other):
nx(other.nx),
nu(other.nu),
ny(other.ny),
A(other.A),
B(other.B),
C(other.C),
D(other.D),
L(other.L),
x(other.x)
{}

// Copy
void zStateObserver::operator=(const zStateObserver& other) {
    if (this != & other) {
        this->nx = other.nx;
        this->nu = other.nu;
        this->ny = other.ny;
        
        this->A = other.A;
        this->B = other.B;
        this->C = other.C;
        this->D = other.D;
        this->L = other.L;

        this->x = other.x;
    }
}

// Move
zStateObserver::zStateObserver(zStateObserver&& other):
nx(std::move(other.nx)),
nu(std::move(other.nu)),
ny(std::move(other.ny)),
A(std::move(other.A)),
B(std::move(other.B)),
C(std::move(other.C)),
D(std::move(other.D)),
L(std::move(other.L)),
x(std::move(other.x))
{}

// Move
void zStateObserver::operator=(zStateObserver&& other) {
    if (this != & other) {
        this->nx = std::move(other.nx);
        this->nu = std::move(other.nu);
        this->ny = std::move(other.ny);
        
        this->A = std::move(other.A);
        this->B = std::move(other.B);
        this->C = std::move(other.C);
        this->D = std::move(other.D);
        this->L = std::move(other.L);

        this->x = std::move(other.x);
    }
}

// Destructor
zStateObserver::~zStateObserver() {
    // Nothing to doo
}



// Getter
// State space dimensions
size_t zStateObserver::get_nx() const {
    return this->nx;
}
size_t zStateObserver::get_nu() const {
    return this->nu;
}
size_t zStateObserver::get_ny() const {
    return this->ny;
}
const Math::Matrix& zStateObserver::get_A() const {
    return this->A;
}
const Math::Matrix& zStateObserver::get_B() const {
    return this->B;
}
const Math::Matrix& zStateObserver::get_C() const {
    return this->C;
}
const Math::Matrix& zStateObserver::get_D() const {
    return this->D;
}
const Math::Matrix& zStateObserver::get_L() const {
    return this->L;
}
const Math::Vector& zStateObserver::get_x() const {
    return this->x;
}


// Setter
void zStateObserver::set_A() {
    this->A = A;
}
void zStateObserver::set_B() {
    this->B = B;
}
void zStateObserver::set_C() {
    this->C = C;
}
void zStateObserver::set_D() {
    this->D = D;
}
void zStateObserver::set_L() {
    this->L = L;
}
void zStateObserver::set_state(const Math::Vector& x0) {
    this->x = x0;
}

// Reset initial state
void zStateObserver::reset() {
    this->x = Math::Vector(this->nx, nullptr);
}




// Discrete-time Luenberger observer
// y^[n] = C * x^[n] + D * u[n]
// e[k] = y[k] - y^[k]
// x^[n+1] = A * x[n] + B * u[n] + L * e[k]
void zStateObserver::update_state(const Math::Vector& u, const Math::Vector& y) {
    const auto e = (y - this->estimated_output(u));
    this->x = this->A * this->x + this->B * u + this->L * e;
}
Math::Vector zStateObserver::estimated_output(const Math::Vector& u) const {
    return this->C * this->x + this->D * u;
}
// calc new output and update internal state
Math::Vector zStateObserver::update(const Math::Vector& u, const Math::Vector& y) {
    this->update_state(u, y);
    return this->state();
}
Math::Vector zStateObserver::operator()(const Math::Vector& u, const Math::Vector& y) {
    return this->update(u, y);
}
Math::Vector zStateObserver::state() const {
    return this->x;
}