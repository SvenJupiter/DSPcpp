#ifndef SJ_DSP_Z_STATE_SPACE_HPP
#define SJ_DSP_Z_STATE_SPACE_HPP

#include "DSP/dsp_types.h"
#include "DSP/Math/Matrix.hpp"
#include "DSP/Math/Vector.hpp"
#include "DSP/Math/Polynomial.hpp"

namespace DSP {

namespace Discrete {


DSP_CLASS zStateSpace {

    protected:
        // State space dimensions
        size_t nx; // Number of states
        size_t nu; // Number of inputs
        size_t ny; // Number of outputs

        Math::Matrix A; // System matrix
        Math::Matrix B; // Input matrix
        Math::Matrix C; // Output matrix
        Math::Matrix D; // Feedthrough matrix
        Math::Vector x; // State vector

        // x[k+1] = A * x[k] + B * u[k]
        // y[k] = C * x[k] + D * u[k]
        void update_state(const Math::Vector& u);
        Math::Vector get_output(const Math::Vector& u) const;

    public:

        // Constructor
        zStateSpace();
        zStateSpace(const size_t num_states, const size_t num_inputs, const size_t num_outputs);
        zStateSpace(const Math::Matrix& A, const Math::Matrix& B);
        zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C);
        zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& D);
        zStateSpace(const Math::Polynomial num, const Math::Polynomial den);

        zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Vector& x0);
        zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Vector& x0);
        zStateSpace(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& D, const Math::Vector& x0);
        zStateSpace(const Math::Polynomial num, const Math::Polynomial den, const Math::Vector& x0);

        // copy
        zStateSpace(const zStateSpace& other);
        void operator=(const zStateSpace& other);

        // move
        zStateSpace(zStateSpace&& other);
        void operator=(zStateSpace&& other);

        // Destructor
        ~zStateSpace();


        // Getter
        // State space dimensions
        size_t get_nx() const;
        size_t get_nu() const;
        size_t get_ny() const;
        const Math::Matrix& get_A() const;
        const Math::Matrix& get_B() const;
        const Math::Matrix& get_C() const;
        const Math::Matrix& get_D() const;
        const Math::Vector& get_x() const;

        // Setter
        void set_A();
        void set_B();
        void set_C();
        void set_D();
        void set_state(const Math::Vector& x0);

        // Reset initial state
        void reset();

        // calc new output and update internal state
        real_t update(const real_t u);
        Math::Vector update(const Math::Vector& u);
        real_t operator()(const real_t& u) { return this->update(u); }
        Math::Vector operator()(const Math::Vector& u) { return this->update(u); }
        Math::Vector output(const Math::Vector& u) const { return this->get_output(u); }
        Math::Vector operator[](const Math::Vector& u) { return this->output(u); }


        enum class sApproximation {
            ForwardEuler,  // s = (z - 1) / Ts
            BackwardEuler, // s = (z - 1) / (z * Ts)
            Trapezoidal // s = (2 / Ts) * (z - 1) / (Z + 1)
        };

        // Static
        static zStateSpace tf2ss(const Math::Polynomial num, const Math::Polynomial den);
        static zStateSpace tf2ss(const Math::Polynomial num, const Math::Polynomial den, const Math::Vector& x0);

        /**
         * @brief Create a PT1 system
         * 
         * @details F(s) = K / (T * s + 1)
         *          F(z) = b / (z - a)
         *          with a = exp(-Ts / T) and b = K * (1 - a)
         * 
         * @param K Gain
         * 
         * @param T Time constant
         * 
         * @param Ts Sample time 
         * 
         * @param x0 Initial state
         * 
         * @return zStateSpace System
         */
        static zStateSpace PT1(const real_t K, const real_t T, const real_t Ts, const real_t x0 = 0);
        
        /**
         * @brief Create a discrete-time lead-lag compensator
         * 
         * @details Continuos: F(s) = (T1 * s + 1) / (T2 * s + 1)
         *          Discrete:  F(z) = (T1 + (Ts - T1) * z^-1) / (T2 + (Ts - T2) * z^-1)
         * 
         * @param T1 Lead time constant
         * 
         * @param T2 Lag time constant
         * 
         * @param Ts Sample time 
         * 
         * @param x0 Initial state
         * 
         * @return zStateSpace System
         */
        static zStateSpace LeadLagCompensator(const real_t T1, const real_t T2, const real_t Ts, const real_t x0 = 0);
        
        /**
         * @brief Create a discrete-time low-pass filter
         * 
         * @details Continuous: F(s) = K / (T * s + 1)
         *          Discrete:   F(z) = K * (Ts / T) * z^-1 / (1 + (Ts / T - 1) * z^-1)
         * 
         * @param K Filter gain
         * 
         * @param T Filter time constant (T = 1 / fc)
         * 
         * @param Ts Sample time 
         * 
         * @param x0 Initial state
         * 
         * @return zStateSpace System
         */
        static zStateSpace LowpassFilter(const real_t K, const real_t T, const real_t Ts, const real_t x0 = 0);
        
        /**
         * @brief Create a Discrete-time washout or high-pass filter
         * 
         * @details Continuous: F(s) = K * T * s / (T * s + 1)
         *          Discrete:   F(z) = K * (1 - z^-1) / (1 + (Ts / T - 1) * z^-1)
         * 
         * @param K Filter gain
         * 
         * @param T Filter time constant (T = 1 / fc)
         * 
         * @param Ts Sample time 
         * 
         * @param x0 Initial state
         * 
         * @return zStateSpace System
         */
        static zStateSpace HighpassFilter(const real_t K, const real_t T, const real_t Ts, const real_t x0 = 0);

        /**
         * @brief Create a Discrete time integrator
         * 
         * @details Continuous: F(s) = 1 / s
         *          Discrete:   
         *              Forward Euler:  F(z) = K * Ts * z^-1 / (1 - z^1)
         *              Backward Euler: F(z) = K * Ts / (1 - z^1)
         *              Trapezoidal:    F(z) = K * (Ts *(1 + z^-1)) / (2 * (1 - z^1))
         * 
         * @param K Integrator gain
         * 
         * @param Ts Sample time 
         * 
         * @param s_approx Approximation: s = f(z)
         * 
         * @param x0 Initial state
         * 
         * @return zStateSpace System
         */
        static zStateSpace Integrator(const real_t K, const real_t Ts, const sApproximation s_approx = sApproximation::ForwardEuler, const real_t x0 = 0);

        /**
         * @brief Create a Discrete time filtered derivative
         * 
         * @details Continuous: F(s) = s / (T * s + 1)
         *          Discrete:   
         *              Forward Euler:  F(z) =  K/T * (1 - z^-1) / (1 + (Ts / T - 1) * z^-1)
         *              Backward Euler: F(z) = K/T * (1 - z^-1) / ((1 + Ts / T) - z^-1)
         *              Trapezoidal:    F(z) = 2*K/T * (1 - z^-1) / ((Ts / T + 2) + (Ts / T - 2) * z^-1)
         * 
         * @param K Derivative gain
         * 
         * @param T Filter time constant (T = 1 / N)
         * 
         * @param Ts Sample time 
         * 
         * @param s_approx Approximation: s = f(z)
         * 
         * @param x0 Initial state
         * 
         * @return zStateSpace System
         */
        static zStateSpace Derivative(const real_t K, const real_t T, const real_t Ts, const sApproximation s_approx = sApproximation::ForwardEuler, const real_t x0 = 0);


}; // class zStateSpace

}; // Namespace Discrete

}; // Namespace DSP

#endif // SJ_DSP_Z_STATE_SPACE_HPP