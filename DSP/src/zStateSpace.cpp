#include <cmath> // expf
#include "DSP/Discrete/zStateSpace.hpp"


using namespace DSP;
using namespace DSP::Discrete;


// Constructor
zStateSpace::zStateSpace():
nx(0), nu(0), ny(0), A(), B(), C(), D(), x()
{}

zStateSpace::zStateSpace(const size_t num_states, const size_t num_inputs, const size_t num_outputs):
nx(num_states), nu(num_inputs), ny(num_outputs),
A(num_states, num_states),
B(num_states, num_inputs),
C(num_outputs, num_states),
D(num_outputs, num_inputs),
x(num_states)
{}

zStateSpace::zStateSpace(const Math::Matrix& A, const Math::Matrix& B):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(A.get_num_rows()), 
A(A), B(B), C(Math::Matrix::eye(A.get_num_rows())), D(Math::Matrix::zeros(A.get_num_rows(), B.get_num_columns())), 
x(A.get_num_rows(), nullptr)
{}

zStateSpace::zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(C.get_num_rows()), 
A(A), B(B), C(C), D(Math::Matrix::zeros(C.get_num_rows(), B.get_num_columns())), 
x(A.get_num_rows(), nullptr)
{}

zStateSpace::zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& D):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(C.get_num_rows()), 
A(A), B(B), C(C), D(D),
x(A.get_num_rows(), nullptr)
{}

zStateSpace::zStateSpace(const Math::Polynomial num, const Math::Polynomial den):
nx(0), nu(0), ny(0), A(), B(), C(), D(), x()
{
    if (num.get_order() == 0 || num.get_order() != den.get_order()) { return; }

    // Init    
    this->nx = den.get_order();
    this->nu = 1;
    this->ny = 1;

    // Create
    this->A = Math::Matrix::zeros(this->nx, this->nx);
    this->B = Math::Matrix::zeros(this->nx, this->nu);
    this->C = Math::Matrix::zeros(this->ny, this->nx);
    this->D = Math::Matrix::zeros(this->ny, this->nu);
    this->x = Math::Vector(this->nx, nullptr);

    // Initilize
    B.m(0, 0) = 1;
    D.m(0, 0) = num.a(0) / den.a(0);
    for (size_t k = 0; k < den.get_order(); ++k) {
        // Init A
        A.m(0, k) = (-1) * den.a(k+1) / den.a(0);
        if (k < den.get_order() - 1) { A.m(k+1, k) = 1; }

        // Init C
        C.m(0, k) = num.a(k+1) / den.a(0) + D.m(0, 0) * A.m(0, k);
    }
}



zStateSpace::zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Vector& x0):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(A.get_num_rows()), 
A(A), B(B), C(Math::Matrix::eye(A.get_num_rows())), D(Math::Matrix::zeros(A.get_num_rows(), B.get_num_columns())), 
x(x0)
{}

zStateSpace::zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Vector& x0):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(C.get_num_rows()), 
A(A), B(B), C(C), D(Math::Matrix::zeros(C.get_num_rows(), B.get_num_columns())), 
x(x0)
{}

zStateSpace::zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& D, const Math::Vector& x0):
nx(A.get_num_rows()), nu(B.get_num_columns()), ny(C.get_num_rows()), 
A(A), B(B), C(C), D(D),
x(x0)
{}

zStateSpace::zStateSpace(const Math::Polynomial num, const Math::Polynomial den, const Math::Vector& x0):
nx(0), nu(0), ny(0), A(), B(), C(), D(), x()
{
    if (num.get_order() == 0 || num.get_order() != den.get_order() || den.get_order() != x0.size()) { return; }

    // Init    
    this->nx = den.get_order();
    this->nu = 1;
    this->ny = 1;

    // Create
    this->A = Math::Matrix::zeros(this->nx, this->nx);
    this->B = Math::Matrix::zeros(this->nx, this->nu);
    this->C = Math::Matrix::zeros(this->ny, this->nx);
    this->D = Math::Matrix::zeros(this->ny, this->nu);
    this->x = x0;

    // Initilize
    B.m(0, 0) = 1;
    D.m(0, 0) = num.a(0) / den.a(0);
    for (size_t k = 0; k < den.get_order(); ++k) {
        // Init A
        A.m(0, k) = (-1) * den.a(k+1) / den.a(0);
        if (k < den.get_order() - 1) { A.m(k+1, k) = 1; }

        // Init C
        C.m(0, k) = num.a(k+1) / den.a(0) + D.m(0, 0) * A.m(0, k);
    }
}




// copy
zStateSpace::zStateSpace(const zStateSpace& other):
nx(other.nx), nu(other.nu), ny(other.ny),
A(other.A), B(other.B), C(other.C), D(other.D),
x(other.x)
{}
void zStateSpace::operator=(const zStateSpace& other) {
    if (this != &other) {
        this->nx = other.nx;
        this->nu = other.nu;
        this->ny = other.ny;
        this->A = other.A;
        this->B = other.B;
        this->C = other.C;
        this->D = other.D;
        this->x = other.x;
    }
}

// move
zStateSpace::zStateSpace(zStateSpace&& other):
nx(std::move(other.nx)), nu(std::move(other.nu)), ny(std::move(other.ny)),
A(std::move(other.A)), B(std::move(other.B)), C(std::move(other.C)), D(std::move(other.D)),
x(std::move(other.x))
{}
void zStateSpace::operator=(zStateSpace&& other) {
    if (this != &other) {
        this->nx = std::move(other.nx);
        this->nu = std::move(other.nu);
        this->ny = std::move(other.ny);
        this->A = std::move(other.A);
        this->B = std::move(other.B);
        this->C = std::move(other.C);
        this->D = std::move(other.D);
        this->x = std::move(other.x);
    }
}

// Destructor
zStateSpace::~zStateSpace() {
    // Nothing to do
}




// Getter
// State space dimensions
size_t zStateSpace::get_nx() const {
    return this->nx;
}
size_t zStateSpace::get_nu() const {
    return this->nu;
}
size_t zStateSpace::get_ny() const {
    return this->ny;
}
const Math::Matrix& zStateSpace::get_A() const {
    return this->A;
}
const Math::Matrix& zStateSpace::get_B() const {
    return this->B;
}
const Math::Matrix& zStateSpace::get_C() const {
    return this->C;
}
const Math::Matrix& zStateSpace::get_D() const {
    return this->D;
}
const Math::Vector& zStateSpace::get_x() const {
    return this->x;
}

// Setter
void zStateSpace::set_A() {

}
void zStateSpace::set_B() {

}
void zStateSpace::set_C() {

}
void zStateSpace::set_D() {

}
void zStateSpace::set_state(const Math::Vector& x0) {
    this->x = x0;
}

// Reset initial state
void zStateSpace::reset() {
    this->x = Math::Vector(this->nx, nullptr);
}

// calc new output and update internal state
real_t zStateSpace::update(const real_t u) {
    return this->update(Math::Vector({u, }))[0];
}
Math::Vector zStateSpace::update(const Math::Vector& u) {
    const auto y = this->get_output(u);
    this->update_state(u);
    return y;
}
Math::Vector zStateSpace::get_output(const Math::Vector& u) const {
    // y[k] = C * x[k] + D * u[k]
    return (this->C * this->x + this->D * u);
}
void zStateSpace::update_state(const Math::Vector& u) {
    // x[k+1] = A * x[k] + B * u[k]
    this->x = this->A * this->x + this->B * u;
}



// Static
zStateSpace zStateSpace::tf2ss(const Math::Polynomial num, const Math::Polynomial den) {
    return zStateSpace(num, den);
}
zStateSpace zStateSpace::tf2ss(const Math::Polynomial num, const Math::Polynomial den, const Math::Vector& x0) {
    return zStateSpace(num, den, x0);
}
zStateSpace zStateSpace::PT1(const real_t K, const real_t T, const real_t Ts, const real_t x0) {

    const real_t a = expf(-Ts / T);
    const real_t b = K * (1 - a);

    const Math::Polynomial num = {0, b};
    const Math::Polynomial den = {1, -a};
    const Math::Vector X0 = {x0, };
    return zStateSpace::tf2ss(num, den, X0);
}
zStateSpace zStateSpace::LeadLagCompensator(const real_t T1, const real_t T2, const real_t Ts, const real_t x0) {

    const Math::Polynomial num = {T1, (Ts - T1)};
    const Math::Polynomial den = {T2, (Ts - T2)};
    const Math::Vector X0 = {x0, };
    return zStateSpace::tf2ss(num, den, X0);
}
zStateSpace zStateSpace::LowpassFilter(const real_t K, const real_t T, const real_t Ts, const real_t x0) {

    const Math::Polynomial num = {0, (K * Ts) / T};
    const Math::Polynomial den = {1, (Ts / T - 1)};
    const Math::Vector X0 = {x0, };
    return zStateSpace::tf2ss(num, den, X0);
}
zStateSpace zStateSpace::HighpassFilter(const real_t K, const real_t T, const real_t Ts, const real_t x0) {

    const Math::Polynomial num = {K, -K};
    const Math::Polynomial den = {1, (Ts / T - 1)};
    const Math::Vector X0 = {x0, };
    return zStateSpace::tf2ss(num, den, X0);
}
zStateSpace zStateSpace::Integrator(const real_t K, const real_t Ts, const sApproximation s_approx, const real_t x0) {

    const Math::Vector X0 = {x0, };
    switch (s_approx) {
        case sApproximation::ForwardEuler: {
            const Math::Polynomial num = {0, K * Ts};
            const Math::Polynomial den = {1, -1};
            return zStateSpace::tf2ss(num, den, X0);
        }
        case sApproximation::BackwardEuler: {
            const Math::Polynomial num = {K * Ts, 0};
            const Math::Polynomial den = {1, -1};
            return zStateSpace::tf2ss(num, den, X0);  
        }
        case sApproximation::Trapezoidal: {
            const Math::Polynomial num = {K * Ts, K * Ts};
            const Math::Polynomial den = {2, -2};
            return zStateSpace::tf2ss(num, den, X0);
        }
        default: {
            const Math::Polynomial num = {0, K * Ts};
            const Math::Polynomial den = {1, -1};
            return zStateSpace::tf2ss(num, den, X0);
        }
    }
}
zStateSpace zStateSpace::Derivative(const real_t K, const real_t T, const real_t Ts, const sApproximation s_approx, const real_t x0) {

    const Math::Vector X0 = {x0, };
    switch (s_approx) {
        case sApproximation::ForwardEuler: {
            const Math::Polynomial num = {K, -K};
            const Math::Polynomial den = {T, (Ts - T)};
            return zStateSpace::tf2ss(num, den, X0);
        }
        case sApproximation::BackwardEuler: {
            const Math::Polynomial num = {K, -K};
            const Math::Polynomial den = {(T+Ts), -T};
            return zStateSpace::tf2ss(num, den, X0);  
        }
        case sApproximation::Trapezoidal: {
            const Math::Polynomial num = {(2 * K), -(2 * K)};
            const Math::Polynomial den = {(Ts + 2 * T), (Ts - 2 * T)};
            return zStateSpace::tf2ss(num, den, X0);
        }
        default: {
            const Math::Polynomial num = {K, -K};
            const Math::Polynomial den = {T, (Ts - T)};
            return zStateSpace::tf2ss(num, den, X0);
        }
    }
}

