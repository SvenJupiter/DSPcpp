#include "DSP/Discrete/pidController.hpp"

using namespace DSP::Discrete;


// Constructor
pidController::pidController():
Ts(0),
Kp(0),
Ki(0),
Kd(0),
Tf(0.001f),
limit_output(false),
upper_limit(1),
lower_limit(0),
anti_windup_method(AntiWindupMethod::None),
Kb(0),
tracking_enabled(false),
Kt(0),
integrator(),
derivative(),
// fromAW(0),
// fromTR(0),
// toInt(0),
preSat(0),
postSat(0)
{}

pidController::pidController(
    const real_t Kp, const real_t Ki, const real_t Kd, const real_t Tf, const real_t Ts, 
    const zStateSpace::sApproximation IF, const zStateSpace::sApproximation DF):
Ts(Ts),
Kp(Kp),
Ki(Ki),
Kd(Kd),
Tf(Tf),
limit_output(false),
upper_limit(1),
lower_limit(0),
anti_windup_method(AntiWindupMethod::None),
Kb(0),
tracking_enabled(false),
Kt(0),
integrator(zStateSpace::Integrator(1, Ts, IF)),
derivative(zStateSpace::Derivative(1, Tf, Ts, DF)),
// fromAW(0),
// fromTR(0),
// toInt(0),
preSat(0),
postSat(0)
{}

// Copy
pidController::pidController(const pidController& other):
Ts(other.Ts),
Kp(other.Kp),
Ki(other.Ki),
Kd(other.Kd),
Tf(other.Tf),
limit_output(other.limit_output),
upper_limit(other.upper_limit),
lower_limit(other.lower_limit),
anti_windup_method(other.anti_windup_method),
Kb(other.Kb),
tracking_enabled(other.tracking_enabled),
Kt(other.Kt),
integrator(other.integrator),
derivative(other.derivative),
// fromAW(other.fromAW),
// fromTR(other.fromTR),
// toInt(other.toInt),
preSat(other.preSat),
postSat(other.postSat)
{}

void pidController::operator=(const pidController& other) {
    if (this != &other) {
        this->Ts = other.Ts;
        this->Kp = other.Kp;
        this->Ki = other.Ki;
        this->Kd = other.Kd;
        this->Tf = other.Tf;
        this->limit_output = other.limit_output;
        this->upper_limit = other.upper_limit;
        this->lower_limit = other.lower_limit;
        this->anti_windup_method = other.anti_windup_method;
        this->Kb = other.Kb;
        this->tracking_enabled = other.tracking_enabled;
        this->Kt = other.Kt;
        this->integrator = other.integrator;
        this->derivative = other.derivative;
        // this->fromAW = other.fromAW;
        // this->fromTR = other.fromTR;
        // this->toInt = other.toInt;
        this->preSat = other.preSat;
        this->postSat = other.postSat;
    }
}

// Move
pidController::pidController(pidController&& other):
Ts(std::move(other.Ts)),
Kp(std::move(other.Kp)),
Ki(std::move(other.Ki)),
Kd(std::move(other.Kd)),
Tf(std::move(other.Tf)),
limit_output(std::move(other.limit_output)),
upper_limit(std::move(other.upper_limit)),
lower_limit(std::move(other.lower_limit)),
anti_windup_method(std::move(other.anti_windup_method)),
Kb(std::move(other.Kb)),
tracking_enabled(std::move(other.tracking_enabled)),
Kt(std::move(other.Kt)),
integrator(std::move(other.integrator)),
derivative(std::move(other.derivative)),
// fromAW(std::move(other.fromAW)),
// fromTR(std::move(other.fromTR)),
// toInt(std::move(other.toInt)),
preSat(std::move(other.preSat)),
postSat(std::move(other.postSat))
{}


void pidController::operator=(pidController&& other) {
    if (this != &other) {
        this->Ts = std::move(other.Ts);
        this->Kp = std::move(other.Kp);
        this->Ki = std::move(other.Ki);
        this->Kd = std::move(other.Kd);
        this->Tf = std::move(other.Tf);
        this->limit_output = std::move(other.limit_output);
        this->upper_limit = std::move(other.upper_limit);
        this->lower_limit = std::move(other.lower_limit);
        this->anti_windup_method = std::move(other.anti_windup_method);
        this->Kb = std::move(other.Kb);
        this->tracking_enabled = std::move(other.tracking_enabled);
        this->Kt = std::move(other.Kt);
        this->integrator = std::move(other.integrator);
        this->derivative = std::move(other.derivative);
        // this->fromAW = std::move(other.fromAW);
        // this->fromTR = std::move(other.fromTR);
        // this->toInt = std::move(other.toInt);
        this->preSat = std::move(other.preSat);
        this->postSat = std::move(other.postSat);
    }
}

// Destructor
pidController::~pidController() {
    // Nothing to do
}


// Setter
void pidController::set_P_gain(const real_t Kp) {
    this->Kp = Kp;
}
void pidController::set_I_gain(const real_t Ki) {
    this->Ki = Ki;
}
void pidController::set_D_gain(const real_t Kd) {
    this->Kd = Kd;
}


// Output saturation
void pidController::set_output_saturation(const bool limit_output, const real_t upper_limit, const real_t lower_limit) {
    this->limit_output = limit_output;
    this->upper_limit = upper_limit;
    this->lower_limit = lower_limit;
}

// Anti-Windup
void pidController::set_anti_windup_method(const AntiWindupMethod method, const real_t Kb) {
    this->anti_windup_method = method;
    this->Kb = Kb;
}

// Enable tracking
void pidController::set_tracking_mode(const bool enable_tracking, const real_t Kt) {
    this->tracking_enabled = enable_tracking;
    this->Kt = Kt;
}


// Calculate new output
real_t pidController::update(const real_t e, const real_t tr) {

    // P
    const real_t fromP = this->Kp * e;

    // I
    const real_t fromTR = (this->tracking_enabled ? this->Kt * (tr - this->postSat) : 0);
    const real_t fromAW = (this->limit_output && this->anti_windup_method == AntiWindupMethod::BackCalculation ? this->Kb * (this->postSat - this->preSat) : 0);
    const real_t preInt = this->Ki * e + fromTR + fromAW;
    const real_t toInt = (
        this->anti_windup_method == AntiWindupMethod::Clamping ? 
        (((this->preSat > this->postSat && preInt > 0) || (this->preSat < this->postSat && preInt < 0)) ? 
        0 : preInt) : preInt);
    const real_t fromI = this->integrator.update(toInt);

    // D
    const real_t fromD = this->derivative.update(this->Kd * e);

    // Saturation
    this->preSat = fromP + fromI + fromD;
    this->postSat = this->preSat;
    if (this->limit_output) {
        if (this->preSat > this->upper_limit) { this->postSat = this->upper_limit; }
        else if (this->preSat < this->lower_limit) { this->postSat = this->lower_limit; }
    }

    // Return
    return this->postSat;
}
real_t pidController::operator()(const real_t e, const real_t tr) {
    return this->update(e, tr);
}

















