#ifndef SJ_DSP_TYPES_H
#define SJ_DSP_TYPES_H

#include <stdint.h> // size_t
#include <stdbool.h> // bool

#ifdef BUILD_SHARDED_DSP_LIB
// #pragma message("Building shared DSP")
#define DSP_FUNCTION __declspec(dllexport)
#endif
#ifdef USE_SHARDED_DSP_LIB
// #pragma message("Using shared DSP")
#define DSP_FUNCTION __declspec(dllimport)
#endif

#ifdef BUILD_STATIC_DSP_LIB
// #pragma message("Building static DSP")
#define DSP_FUNCTION
#endif
#ifdef USE_STATIC_DSP_LIB
// #pragma message("Using static DSP")
#define DSP_FUNCTION
#endif

#ifndef DSP_FUNCTION
#define DSP_FUNCTION
#endif

#define DSP_TYPE DSP_FUNCTION
#define DSP_CLASS class DSP_TYPE 
#define DSP_ENUM_CLASS enum class DSP_TYPE 

typedef float real_t;



#endif // SJ_DSP_TYPES_H