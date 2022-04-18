#include <limits>
#include <cstdint>
#include <utility>
#include <cassert>

namespace pbrt{
namespace gpu{

class Scene;
class GeometricPrimitive;
class SurfaceInteraction;
class Shape;
class Sphere;
class Medium;
struct MediumInterface;
template <int nSpectrumSamples>
class CoefficientSpectrum;
class RGBSpectrum;
class SampledSpectrum;
#ifdef PBRT_SAMPLED_SPECTRUM
  typedef SampledSpectrum Spectrum;
#else
  typedef RGBSpectrum Spectrum;
#endif

template <typename T>
class Vector2;
template <typename T>
class Vector3;
template <typename T>
class Point3;
template <typename T>
class Point2;
template <typename T>
class Normal3;
class Ray;
class RayDifferential;
template <typename T>
class Bounds2;
template <typename T>
class Bounds3;

class Sampler;

#ifdef PBRT_FLOAT_AS_DOUBLE
  typedef double Float;
#else
  typedef float Float;
#endif  // PBRT_FLOAT_AS_DOUBLE

// Global Constants
// extern PBRT_CONSTEXPR Float MaxFloat = std::numeric_limits<Float>::max();
// PBRT_CONSTEXPR Float Infinity = std::numeric_limits<Float>::infinity();
// PBRT_CONSTEXPR Float MachineEpsilon =
//     std::numeric_limits<Float>::epsilon() * 0.5;
// PBRT_CONSTEXPR Float ShadowEpsilon = 0.0001f;
// PBRT_CONSTEXPR Float Pi = 3.14159265358979323846;
// PBRT_CONSTEXPR Float InvPi = 0.31830988618379067154;
// PBRT_CONSTEXPR Float Inv2Pi = 0.15915494309189533577;
// PBRT_CONSTEXPR Float Inv4Pi = 0.07957747154594766788;
// PBRT_CONSTEXPR Float PiOver2 = 1.57079632679489661923;
// PBRT_CONSTEXPR Float PiOver4 = 0.78539816339744830961;
// PBRT_CONSTEXPR Float Sqrt2 = 1.41421356237309504880;

// class Scene;
// class Integrator;
// class SamplerIntegrator;

// class Transform;
// struct Interaction;
// class SurfaceInteraction;
// class Shape;
// class Primitive;
// class GeometricPrimitive;
// class TransformedPrimitive;
// class Camera;
// struct CameraSample;
// class ProjectiveCamera;
// class Filter;
// class Film;
// class FilmTile;
// class BxDF;
// class BRDF;
// class BTDF;
// class BSDF;
// class Material;
// template <typename T>
// class Texture;
// class Medium;
// class MediumInteraction;
// struct MediumInterface;
// class BSSRDF;
// class SeparableBSSRDF;
// class TabulatedBSSRDF;
// struct BSSRDFTable;
// class Light;
// class VisibilityTester;
// class AreaLight;
// struct Distribution1D;
// class Distribution2D;
// #ifdef PBRT_FLOAT_AS_DOUBLE
//   typedef double Float;
// #else
//   typedef float Float;
// #endif  // PBRT_FLOAT_AS_DOUBLE
// class RNG;
// class ProgressReporter;
// class MemoryArena;
// template <typename T, int logBlockSize = 2>
// class BlockedArray;
// struct Matrix4x4;
// class ParamSet;
// template <typename T>
// struct ParamSetItem;
// struct Options {
//     Options() {
//         cropWindow[0][0] = 0;
//         cropWindow[0][1] = 1;
//         cropWindow[1][0] = 0;
//         cropWindow[1][1] = 1;
//     }
//     int nThreads = 0;
//     bool quickRender = false;
//     bool quiet = false;
//     bool cat = false, toPly = false;
//     std::string imageFile;
//     // x0, x1, y0, y1
//     Float cropWindow[2][2];
// };

// extern Options PbrtOptions;
// class TextureParams;

__device__
Spectrum Li(const RayDifferential &r, const Scene &scene, Sampler &sampler, int depth);

// Global Inline Functions
__device__
uint32_t FloatToBits(float f) {
    uint32_t ui;
    memcpy(&ui, &f, sizeof(float));
    return ui;
}

__device__
float BitsToFloat(uint32_t ui) {
    float f;
    memcpy(&f, &ui, sizeof(uint32_t));
    return f;
}

__device__
uint64_t FloatToBits(double f) {
    uint64_t ui;
    memcpy(&ui, &f, sizeof(double));
    return ui;
}

__device__
double BitsToFloat(uint64_t ui) {
    double f;
    memcpy(&f, &ui, sizeof(uint64_t));
    return f;
}

__device__
float NextFloatUp(float v) {
    // Handle infinity and negative zero for _NextFloatUp()_
    if (isinf(v) && v > 0.) return v;
    if (v == -0.f) v = 0.f;

    // Advance _v_ to next higher float
    uint32_t ui = FloatToBits(v);
    if (v >= 0)
        ++ui;
    else
        --ui;
    return BitsToFloat(ui);
}

__device__
float NextFloatDown(float v) {
    // Handle infinity and positive zero for _NextFloatDown()_
    if (isinf(v) && v < 0.) return v;
    if (v == 0.f) v = -0.f;
    uint32_t ui = FloatToBits(v);
    if (v > 0)
        --ui;
    else
        ++ui;
    return BitsToFloat(ui);
}

__device__
double NextFloatUp(double v, int delta = 1) {
    if (isinf(v) && v > 0.) return v;
    if (v == -0.f) v = 0.f;
    uint64_t ui = FloatToBits(v);
    if (v >= 0.)
        ui += delta;
    else
        ui -= delta;
    return BitsToFloat(ui);
}

__device__
double NextFloatDown(double v, int delta = 1) {
    if (isinf(v) && v < 0.) return v;
    if (v == 0.f) v = -0.f;
    uint64_t ui = FloatToBits(v);
    if (v > 0.)
        ui -= delta;
    else
        ui += delta;
    return BitsToFloat(ui);
}

__device__
Float gamma(int n) {
    return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

__device__
Float GammaCorrect(Float value) {
    if (value <= 0.0031308f) return 12.92f * value;
    return 1.055f * std::pow(value, (Float)(1.f / 2.4f)) - 0.055f;
}

__device__
Float InverseGammaCorrect(Float value) {
    if (value <= 0.04045f) return value * 1.f / 12.92f;
    return std::pow((value + 0.055f) * 1.f / 1.055f, (Float)2.4f);
}

template <typename T, typename U, typename V>
__device__
T Clamp(T val, U low, V high) {
    if (val < low)
        return low;
    else if (val > high)
        return high;
    else
        return val;
}

template <typename T>
__device__
T Mod(T a, T b) {
    T result = a - (a / b) * b;
    return (T)((result < 0) ? result + b : result);
}

template <>
__device__
Float Mod(Float a, Float b) {
    return std::fmod(a, b);
}

__device__
Float Radians(Float deg) { return (Pi / 180) * deg; }

__device__
Float Degrees(Float rad) { return (180 / Pi) * rad; }

__device__
Float Log2(Float x) {
    const Float invLog2 = 1.442695040888963387004650940071;
    return std::log(x) * invLog2;
}

__device__
int Log2Int(uint32_t v) {
    assert(false);
}

__device__
int Log2Int(int32_t v) { return Log2Int((uint32_t)v); }

__device__
int Log2Int(uint64_t v) {
    assert(false);
}

__device__
int Log2Int(int64_t v) { return Log2Int((uint64_t)v); }

template <typename T>
__device__
PBRT_CONSTEXPR bool IsPowerOf2(T v) {
    return v && !(v & (v - 1));
}

__device__
int32_t RoundUpPow2(int32_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return v + 1;
}

__device__
int64_t RoundUpPow2(int64_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    return v + 1;
}

__device__
int CountTrailingZeros(uint32_t v) {
    assert(false);
}

template <typename Predicate>
__device__
int FindInterval(int size, const Predicate &pred) {
    int first = 0, len = size;
    while (len > 0) {
        int half = len >> 1, middle = first + half;
        // Bisect range based on value of _pred_ at _middle_
        if (pred(middle)) {
            first = middle + 1;
            len -= half + 1;
        } else
            len = half;
    }
    return Clamp(first - 1, 0, size - 2);
}

__device__
Float Lerp(Float t, Float v1, Float v2) { return (1 - t) * v1 + t * v2; }

__device__
bool Quadratic(Float a, Float b, Float c, Float *t0, Float *t1) {
    // Find quadratic discriminant
    double discrim = (double)b * (double)b - 4 * (double)a * (double)c;
    if (discrim < 0) return false;
    double rootDiscrim = std::sqrt(discrim);

    // Compute quadratic _t_ values
    double q;
    if (b < 0)
        q = -.5 * (b - rootDiscrim);
    else
        q = -.5 * (b + rootDiscrim);
    *t0 = q / a;
    *t1 = c / q;
    if (*t0 > *t1) {Float tmp = *t0; *t0 = *t1; *t1 = tmp;}
    return true;
}

__device__
Float ErfInv(Float x) {
    Float w, p;
    x = Clamp(x, -.99999f, .99999f);
    w = -std::log((1 - x) * (1 + x));
    if (w < 5) {
        w = w - 2.5f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p * w;
        p = -3.5233877e-06f + p * w;
        p = -4.39150654e-06f + p * w;
        p = 0.00021858087f + p * w;
        p = -0.00125372503f + p * w;
        p = -0.00417768164f + p * w;
        p = 0.246640727f + p * w;
        p = 1.50140941f + p * w;
    } else {
        w = std::sqrt(w) - 3;
        p = -0.000200214257f;
        p = 0.000100950558f + p * w;
        p = 0.00134934322f + p * w;
        p = -0.00367342844f + p * w;
        p = 0.00573950773f + p * w;
        p = -0.0076224613f + p * w;
        p = 0.00943887047f + p * w;
        p = 1.00167406f + p * w;
        p = 2.83297682f + p * w;
    }
    return p * x;
}

__device__
Float Erf(Float x) {
    // constants
    Float a1 = 0.254829592f;
    Float a2 = -0.284496736f;
    Float a3 = 1.421413741f;
    Float a4 = -1.453152027f;
    Float a5 = 1.061405429f;
    Float p = 0.3275911f;

    // Save the sign of x
    int sign = 1;
    if (x < 0) sign = -1;
    x = std::abs(x);

    // A&S formula 7.1.26
    Float t = 1 / (1 + p * x);
    Float y =
        1 -
        (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

    return sign * y;
}

};
};