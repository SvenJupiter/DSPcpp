
# define Library "DSP"
target_sources(DSPcpp PRIVATE 
    Polynomial.cpp 
    Matrix.cpp
    Vector.cpp
    Signal.cpp
    zStateSpace.cpp
    zStateObserver.cpp
    pidController.cpp
)
