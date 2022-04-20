#ifndef SJ_DSP_Z_STATE_OBSERVER_HPP
#define SJ_DSP_Z_STATE_OBSERVER_HPP

#include "DSP/dsp_types.h"
#include "DSP/Discrete/zStateSpace.hpp"

namespace DSP {

namespace Discrete {


// Discrete-time Luenberger observer
DSP_CLASS zStateObserver {

    protected:
        // State space dimensions
        size_t nx; // Number of states
        size_t nu; // Number of inputs
        size_t ny; // Number of outputs

        Math::Matrix A; // System matrix
        Math::Matrix B; // Input matrix
        Math::Matrix C; // Output matrix
        Math::Matrix D; // Feedthrough matrix
        Math::Matrix L; // Observer gain matrix
        Math::Vector x; // Estimated state vector

        // Discrete-time Luenberger observer
        // y^[n] = C * x^[n] + D * u[n]
        // e[k] = y[k] - y^[k]
        // x^[n+1] = A * x^[n] + B * u[n] + L * e[k]
        void update_state(const Math::Vector& u, const Math::Vector& y);
        Math::Vector estimated_output(const Math::Vector& u) const;

    public:

        // Constructor
        zStateObserver();
        zStateObserver(const size_t num_states, const size_t num_inputs, const size_t num_outputs);
        zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& L);
        zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& L);
        zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& D, const Math::Matrix& L);
        zStateObserver(const Math::Polynomial& num, const Math::Polynomial& den, const Math::Matrix& L);
        zStateObserver(const zStateSpace& zss, const Math::Matrix& L);
        zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& L, const Math::Vector& x0);
        zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& L, const Math::Vector& x0);
        zStateObserver(const Math::Matrix& A, const Math::Matrix& B, const Math::Matrix& C, const Math::Matrix& D, const Math::Matrix& L, const Math::Vector& x0);
        zStateObserver(const Math::Polynomial& num, const Math::Polynomial& den, const Math::Matrix& L, const Math::Vector& x0);
        zStateObserver(const zStateSpace& zss, const Math::Matrix& L, const Math::Vector& x0);


        // Copy
        zStateObserver(const zStateObserver& other);
        void operator=(const zStateObserver& other);

        // Move
        zStateObserver(zStateObserver&& other);
        void operator=(zStateObserver&& other);

        // Destructor
        virtual ~zStateObserver();

        // Getter
        // State space dimensions
        size_t get_nx() const;
        size_t get_nu() const;
        size_t get_ny() const;
        const Math::Matrix& get_A() const;
        const Math::Matrix& get_B() const;
        const Math::Matrix& get_C() const;
        const Math::Matrix& get_D() const;
        const Math::Matrix& get_L() const;
        const Math::Vector& get_x() const;

        // Setter
        void set_A();
        void set_B();
        void set_C();
        void set_D();
        void set_L();
        void set_state(const Math::Vector& x0);

        // Reset initial state
        void reset();

        // calc new output and update internal state
        Math::Vector update(const Math::Vector& u, const Math::Vector& y);
        Math::Vector operator()(const Math::Vector& u, const Math::Vector& y);
        Math::Vector state() const;


}; // class zStateObserver

}; // Namespace Discrete

}; // Namespace DSP

#endif // SJ_DSP_Z_STATE_OBSERVER_HPP