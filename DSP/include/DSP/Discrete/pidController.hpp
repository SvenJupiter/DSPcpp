#ifndef SJ_PID_CONTROLLER_H
#define SJ_PID_CONTROLLER_H

#include "DSP/Discrete/zStateSpace.hpp"

namespace DSP {

namespace Discrete {

DSP_CLASS pidController {

    public:
        // Anti-Windup
        enum class AntiWindupMethod {
            None,
            BackCalculation,
            Clamping
        };

    protected:

        // Sample time
        real_t Ts;

        // Controller Gain
        real_t Kp;
        real_t Ki;
        real_t Kd;
        real_t Tf;

        // Output saturation
        bool limit_output;
        real_t upper_limit;
        real_t lower_limit;

        // Anti-Windup
        AntiWindupMethod anti_windup_method;
        real_t Kb;
        
        // Tracking
        bool tracking_enabled;
        real_t Kt;

        // Integrator & Derivative
        zStateSpace integrator;
        zStateSpace derivative;

        // internal states
        // real_t fromAW;
        // real_t fromTR;
        // real_t toInt;
        real_t preSat;
        real_t postSat;
        // real_t output;



    public:
        // Constructor
        pidController();
        pidController(
            const real_t Kp, const real_t Ki, const real_t Kd, const real_t Tf, const real_t Ts, 
            const zStateSpace::sApproximation IF = zStateSpace::sApproximation::ForwardEuler, 
            const zStateSpace::sApproximation DF = zStateSpace::sApproximation::ForwardEuler
        );

        // Copy
        pidController(const pidController& other);
        void operator=(const pidController& other);

        // Move
        pidController(pidController&& other);
        void operator=(pidController&& other);

        // Destructor
        virtual ~pidController();


        // Setter
        void set_P_gain(const real_t Kp);
        void set_I_gain(const real_t Ki);
        void set_D_gain(const real_t Kd);

        // Output saturation
        void set_output_saturation(const bool limit_output, const real_t upper_limit = 1, const real_t lower_limit = 0);

        // Anti-Windup
        void set_anti_windup_method(const AntiWindupMethod method, const real_t Kb = 1);

        // Enable tracking
        void set_tracking_mode(const bool enable_tracking, const real_t Kt = 1);


        // Calculate new output
        real_t update(const real_t e, const real_t tr = 0);
        real_t operator()(const real_t e, const real_t tr = 0);


}; // class pidController

}; // namespace Discrete

}; // namespace DSP


#endif // SJ_PID_CONTROLLER_H